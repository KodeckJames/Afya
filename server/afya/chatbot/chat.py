"""
## Build Medical Question Answering System using LangChain and Gemini Models
"""

# Uncomment and install required libraries if not already done


import os
from langchain.document_loaders import JSONLoader
from langchain.text_splitter import TokenTextSplitter
from langchain_google_genai import GoogleGenerativeAIEmbeddings,ChatGoogleGenerativeAI
from langchain.vectorstores import FAISS
from langchain.chains import RetrievalQA
from langchain.prompts import PromptTemplate
from dotenv import load_dotenv
load_dotenv()

google_api_key=os.getenv("GOOGLE_API_KEY")


# Load PubMed articles from a JSON file
def metadata_func(record: dict, metadata: dict) -> dict:
    metadata["year"] = record.get("pub_date").get('year')
    metadata["month"] = record.get("pub_date").get('month')
    metadata["day"] = record.get("pub_date").get('day')
    metadata["title"] = record.get("article_title")
    return metadata

loader = JSONLoader(
    file_path='output.json',
    jq_schema='.[]',
    content_key='article_abstract',
    metadata_func=metadata_func
)
data = loader.load()
print(f"{len(data)} pubmed articles are loaded!")
print(data[1])

# Chunk abstracts into small text passages for efficient retrieval and LLM context length
text_splitter = TokenTextSplitter(chunk_size=128, chunk_overlap=64)
chunks = text_splitter.split_documents(data)
print(f"{len(data)} pubmed articles are converted to {len(chunks)} text fragments!")
print(chunks[0])

embeddings = GoogleGenerativeAIEmbeddings(
    model="models/text-embedding-004",
    api_key=google_api_key,
)

# Build the vector database (VDB) using FAISS
db = FAISS.from_documents(chunks, embeddings)

llm=ChatGoogleGenerativeAI(
    model="gemini-1.5-flash",
    api_key=google_api_key,
)

# Define the RAG pipeline using LangChain
PROMPT_TEMPLATE_2 = """You are a medical assistant for question-answering tasks. Answer the Question using the provided Context only. Your answer should be in your own words and be no longer than 128 words.

Context: {context}

Question: {question}

Answer:"""
PROMPT2 = PromptTemplate.from_template(PROMPT_TEMPLATE_2)

db_retriever = db.as_retriever(k=2)
qa_chain = RetrievalQA.from_chain_type(
    llm,
    retriever=db_retriever,
    chain_type_kwargs={"prompt": PROMPT2},
    return_source_documents=True
)

# Example usage
while True:
    question = input("Enter your medical question (or type 'exit' to quit): ")
    if question.lower() == 'exit':
        break
    result = qa_chain({"query": question})
    print("\nAnswer:", result['result'])
    print("\nSources:")
    for doc in result['source_documents']:
        print(f"- {doc.metadata['title']} ({doc.metadata['year']})")
