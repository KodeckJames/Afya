## Extract medical data to be used for the chatbot
from typing import List, Dict
from Bio import Entrez
from langchain_community.vectorstores import Chroma
from langchain.text_splitter import RecursiveCharacterTextSplitter
from langchain.prompts import PromptTemplate
from langchain_google_genai import ChatGoogleGenerativeAI,GoogleGenerativeAIEmbeddings
from langchain.chains import RetrievalQA
from tqdm import tqdm
import time

## Create a pubmed RAG Q&A chatbot
class PubmedRAG:
    def __init__(self,email: str,med_api_key: str,google_api_key: str):
        """
                Initialize the PubMed RAG system

                Args:
                    email: Email for PubMed API
                    med_api_key: PubMed API key
                    google_api_key:GEMINI API key for embeddings and generation
                """

        self.email = email
        Entrez.email = email
        self.api_key = med_api_key
        if med_api_key:
            Entrez.api_key = med_api_key

        ## Initialize embedding and LLM
        self.embeddings = GoogleGenerativeAIEmbeddings(
            model="models/text-embedding-004",
            api_key = google_api_key
        )

        self.llm = ChatGoogleGenerativeAI(
            model='gemini-1.5-flash',
            verbose=True,
            google_api_key= google_api_key
        )

        # Initialize text splitter
        self.text_splitter = RecursiveCharacterTextSplitter(
            chunk_size=1000,
            chunk_overlap=200,
            separators=["\n\n", "\n", " ", ""]
        )

        # Initialize vector store
        self.vectorstore = None

    def fetch_pubmed_articles(self, query: str, max_results: int = 100) -> List[Dict]:
        """
        Fetch articles from PubMed based on search query
        """
        # Search for article IDs
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
        record = Entrez.read(handle)
        handle.close()

        # Fetch article details
        articles = []
        for i in tqdm(range(0, len(record["IdList"]), 50)):
            batch = record["IdList"][i:i + 50]
            handle = Entrez.efetch(db="pubmed", id=','.join(batch), rettype="medline", retmode="text")
            articles.extend(self._parse_medline(handle.read()))
            handle.close()
            time.sleep(0.34)  # Rate limiting

        return articles

    def _parse_medline(self, medline_text: str) -> List[Dict]:
        """
         format text into structured data
        """
        articles = []
        current_article = {}

        for line in medline_text.split('\n'):
            if line.startswith('PMID-'):
                if current_article:
                    articles.append(current_article)
                current_article = {'pmid': line[6:].strip()}
            elif line.startswith('TI  -'):
                current_article['title'] = line[6:].strip()
            elif line.startswith('AB  -'):
                current_article['abstract'] = line[6:].strip()
            elif line.startswith('     '):  # Continuation of previous field
                if 'abstract' in current_article:
                    current_article['abstract'] += ' ' + line.strip()
                elif 'title' in current_article:
                    current_article['title'] += ' ' + line.strip()

        if current_article:
            articles.append(current_article)

        return articles

    def create_knowledge_base(self, articles: List[Dict]):
        """
        Create vector store from fetched articles
        """
        # Prepare documents for vectorization
        documents = []
        for article in articles:
            if 'abstract' in article:
                text = f"Title: {article.get('title', '')}\nAbstract: {article['abstract']}\nPMID: {article['pmid']}"
                chunks = self.text_splitter.split_text(text)
                documents.extend(chunks)

        # Create vector store
        self.vectorstore = Chroma.from_texts(
            texts=documents,
            embedding=self.embeddings,
            persist_directory="./pubmed_chroma_db"
        )

    def query_knowledge_base(self, query: str, k: int = 4) -> str:
        """
        Query the knowledge base using RAG
        """
        if not self.vectorstore:
            raise ValueError("Knowledge base not initialized. Call create_knowledge_base first.")

        # Create RAG chain
        prompt_template = """
            You are a friendly medical AI doctor assistant helping to answer questions from patients with various illnesses based on scientific literature.
            Use the following pieces of context from medical research papers to answer the question.
            If you cannot answer the question based on the context, say so.
            Always include PMID references for your sources.

            Context: {context}

            Question: {question}

            Answer: """

        PROMPT = PromptTemplate(
            template=prompt_template,
            input_variables=["context", "question"]
        )

        chain = RetrievalQA.from_chain_type(
            llm=self.llm,
            chain_type="stuff",
            retriever=self.vectorstore.as_retriever(search_kwargs={"k": k}),
            chain_type_kwargs={"prompt": PROMPT}
        )

        return chain.invoke(query)

    def main():
        # Example usage
        email = "your_email@example.com"
        med_api_key = "your_pubmed_api_key"  # Optional
        google_api_key = "your_gemini_api_key"

        # Initialize system
        rag_system = PubmedRAG(email, med_api_key, google_api_key)

        # Fetch articles
        articles = rag_system.fetch_pubmed_articles(
            query="cancer immunotherapy treatment outcomes",
            max_results=100
        )

        # Create knowledge base
        rag_system.create_knowledge_base(articles)

        # Query the system
        query = "What are the latest developments in CAR-T cell therapy for solid tumors?"
        response = rag_system.query_knowledge_base(query)
        print(response)

    if __name__ == "__main__":
        main()
