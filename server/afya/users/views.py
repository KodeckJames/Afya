from django.shortcuts import render
from rest_framework import generics
from .serializers import UserSerializer
from .models import CustomUser

class UserListCreateView(generics.ListCreateAPIView):
    queryset = CustomUser.objects.all()
    serializer_class = UserSerializer
