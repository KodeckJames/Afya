import os
import time
import random
from Bio import Entrez
from langchain_google_genai import GoogleGenerativeAIEmbeddings,ChatGoogleGenerativeAI
from langchain.text_splitter import RecursiveCharacterTextSplitter
from langchain_community.vectorstores import Chroma
from langchain.prompts import PromptTemplate
from dotenv import load_dotenv

load_dotenv()


class PubmedRAG:
    """
    Class to query PubMed, build a knowledge base using LangChain, and provide retrieval-augmented responses.
    """

    def __init__(self, email, med_api_key, google_api_key, collection_name="pubmed_knowledge"):
        # Initialize PubMed API access
        Entrez.email = email
        Entrez.api_key = med_api_key

        # Initialize Google GenAI for embedding and LLM
        self.google_api_key = google_api_key
        self.embeddings = GoogleGenerativeAIEmbeddings(model="models/text-embedding-004",google_api_key=self.google_api_key)
        self.llm = ChatGoogleGenerativeAI(model='gemini-1.5-flash',google_api_key=self.google_api_key)

        # Vector store configuration
        self.collection_name = collection_name
        self.vectorstore = Chroma(
            collection_name=self.collection_name,
            embedding_function=self.embeddings
        )

        # Conversation state
        self.history = []

    def fetch_pubmed_articles(self, query, max_results=10, min_abstract_length=100):
        """Fetch articles from PubMed based on a query."""
        try:
            handle = Entrez.esearch(
                db="pubmed", term=query, retmax=max_results
            )
            record = Entrez.read(handle)
            handle.close()
            pmids = record["IdList"]

            articles = []
            for pmid in pmids:
                handle = Entrez.efetch(
                    db="pubmed", id=pmid, rettype="medline", retmode="text"
                )
                medline_text = handle.read()
                handle.close()

                article = self._parse_medline(medline_text)
                if len(article["abstract"]) >= min_abstract_length:
                    articles.append(article)

            return articles
        except Exception as e:
            print(f"Error fetching PubMed articles: {e}")
            return []

    def _parse_medline(self, medline_text):
        """Parse MEDLINE formatted text to extract article information."""
        article = {"title": "", "abstract": "", "authors": [], "mesh_terms": []}
        lines = medline_text.splitlines()

        for line in lines:
            if line.startswith("TI  -"):
                article["title"] = line[5:].strip()
            elif line.startswith("AB  -"):
                article["abstract"] += line[5:].strip() + " "
            elif line.startswith("AU  -"):
                article["authors"].append(line[5:].strip())
            elif line.startswith("MH  -"):
                article["mesh_terms"].append(line[5:].strip())

        article["abstract"] = article["abstract"].strip()
        return article

    def create_knowledge_base(self, articles, batch_size=5):
        """Build and store embeddings for PubMed articles in batches."""
        try:
            text_splitter = RecursiveCharacterTextSplitter(
                chunk_size=1000, chunk_overlap=100
            )

            for i in range(0, len(articles), batch_size):
                batch = articles[i:i + batch_size]
                docs = []

                for article in batch:
                    content = f"{article['title']}\n{article['abstract']}"
                    splits = text_splitter.split_text(content)
                    for split in splits:
                        docs.append({"content": split, "metadata": article})

                # Add documents to vector store
                retries = 3
                for attempt in range(retries):
                    try:
                        self.vectorstore.add_texts([doc["content"] for doc in docs],
                                                   metadata=[doc["metadata"] for doc in docs])
                        break
                    except Exception as e:
                        if attempt < retries - 1:
                            print(f"Error adding to vector store, retrying ({attempt + 1}/{retries}): {e}")
                            time.sleep(random.uniform(1, 3))
                        else:
                            print(f"Failed to add texts after {retries} attempts: {e}")

        except Exception as e:
            print(f"Error creating knowledge base: {e}")

    def query_knowledge_base(self, query, k=3):
        """Query the vector store and generate a response using LLM."""
        try:
            # Perform vector store similarity search
            results = self.vectorstore.similarity_search(query, k=k)
            context = "\n\n".join([res["content"] for res in results])

            prompt_template = PromptTemplate(
                template="""
                You are an AI assistant using PubMed data to answer questions. Here is some context to assist:

                {context}

                Question: {query}

                Provide a detailed, accurate response based on the context above.
                """
            )
            prompt = prompt_template.format(context=context, query=query)
            response = self.llm(prompt)

            # Record exchange history
            self.history.append({"query": query, "response": response})

            return response
        except Exception as e:
            print(f"Error querying knowledge base: {e}")
            return "An error occurred while processing your request."

    def get_conversation_context(self, max_exchanges=5):
        """Return the most recent exchanges for context."""
        return self.history[-max_exchanges:]

    def clear_history(self):
        """Clear conversation history."""
        self.history = []

    def save_session(self, filename):
        """Save conversation history to a file."""
        try:
            with open(filename, "w") as file:
                for entry in self.history:
                    file.write(f"Query: {entry['query']}\nResponse: {entry['response']}\n\n")
        except Exception as e:
            print(f"Error saving session: {e}")

    def get_session_summary(self):
        """Return a summary of the session metrics."""
        total_exchanges = len(self.history)
        success_rate = sum(1 for h in self.history if h[
            'response'] != "An error occurred while processing your request.") / total_exchanges if total_exchanges > 0 else 0

        return {
            "total_exchanges": total_exchanges,
            "success_rate": success_rate,
            "session_duration": f"{time.time() - self.start_time:.2f} seconds"
        }


def main():
    ## Initialize the PubmedRAG class
    email = "otienookoth007@gmail.com",
    med_api_key = os.getenv("MED_API_KEY"),
    google_api_key = os.getenv("GOOGLE_API_KEY")

    pubmed_rag = PubmedRAG(email,med_api_key,google_api_key)

    articles = pubmed_rag.fetch_pubmed_articles("COVID-19", max_results=5)
    pubmed_rag.create_knowledge_base(articles)

    response = pubmed_rag.query_knowledge_base("What are the recent advancements in COVID-19 vaccines?")
    print(response)

if __name__ == "__main__":
    main()
