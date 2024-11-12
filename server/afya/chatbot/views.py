from django.shortcuts import render
from rest_framework.response import Response
from rest_framework.decorators import api_view
from .pubmed_rag import PubmedRAG  # Adjust the import based on your file structure

@api_view(['POST'])
def chat(request):
    user_message = request.data.get("message")
    if not user_message:
        return Response({"error": "No message provided"}, status=400)

    # Instantiate PubmedRAG (ensure API keys are properly configured)
    pubmed_rag = PubmedRAG(
        email="otienookoth007@gmail.com",
        med_api_key=os.getenv("MED_API_KEY"),
        google_api_key=os.getenv("GOOGLE_API_KEY")
    )

    # Generate a response from PubmedRAG
    response_text = pubmed_rag.query_knowledge_base(user_message)
    
    return Response({"response": response_text})

