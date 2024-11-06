import React, { useState, useRef, useEffect } from 'react'
import {
  View,
  Text,
  TextInput,
  TouchableOpacity,
  FlatList,
  KeyboardAvoidingView,
} from 'react-native'
import { Ionicons } from '@expo/vector-icons' // For the send icon

// Define the Message interface
interface Message {
  text: string
  isUser: boolean
}

const ChatScreen: React.FC = () => {
  // State for messages and input text
  const [messages, setMessages] = useState<Message[]>([])
  const [inputText, setInputText] = useState<string>('')

  // States for the welcome messages
  const [welcomeMessage, setWelcomeMessage] = useState<string>('')
  const [helpMessage, setHelpMessage] = useState<string>('')

  // Create a ref for FlatList to access its scroll functions
  const flatListRef = useRef<FlatList>(null)

  // Scroll to the bottom when the messages state changes
  useEffect(() => {
    flatListRef.current?.scrollToEnd({ animated: true })
  }, [messages])

  // Simulate typing effect for welcome and help messages
  useEffect(() => {
    let messageTimeout: NodeJS.Timeout
    const typeMessage = (
      message: string,
      setter: React.Dispatch<React.SetStateAction<string>>,
      delay: number
    ) => {
      let index = 0
      const interval = setInterval(() => {
        setter((prev) => prev + message[index])
        index += 1
        if (index === message.length) {
          clearInterval(interval)
        }
      }, delay)
    }

    // Start typing "Hi there" and "How can I help you today?" when the component mounts
    typeMessage('Hi there 👋', setWelcomeMessage, 100)
    setTimeout(() => {
      typeMessage('How can I help you today?', setHelpMessage, 100)
    }, 1500)

    return () => clearTimeout(messageTimeout)
  }, [])

  // Stop typing effect once the user starts typing
  const handleUserInput = () => {
    setMessages((prevMessages) => [
      ...prevMessages,
      { text: inputText, isUser: true },
    ])
    setInputText('') // Clear input
    // Add the AI response after a short delay
    setTimeout(() => {
      setMessages((prevMessages) => [
        ...prevMessages,
        { text: 'AI Response', isUser: false },
      ])
    }, 1000)
  }

  return (
    <KeyboardAvoidingView className="flex-1 bg-gray-900" behavior="padding">
      <FlatList
        ref={flatListRef} // Assign the ref to FlatList
        data={messages}
        renderItem={({ item }) => (
          <View
            className={`my-2 mx-4 p-3 rounded-lg max-w-[75%] ${
              item.isUser
                ? 'self-end bg-gradient-to-r from-purple-500 to-indigo-600'
                : 'self-start bg-gray-700'
            }`}
          >
            <Text className={`${item.isUser ? 'text-white' : 'text-gray-200'}`}>
              {item.text}
            </Text>
          </View>
        )}
        keyExtractor={(item, index) => index.toString()}
      />

      {/* Welcome Message Area */}
      {welcomeMessage || helpMessage ? (
        <View className="flex-1 items-center">
          {welcomeMessage ? (
            <Text className="text-purple-300 text-4xl ">{welcomeMessage}</Text>
          ) : null}
          {helpMessage ? (
            <Text className="text-purple-300 text-3xl">{helpMessage}</Text>
          ) : null}
        </View>
      ) : null}

      {/* Input Area */}
      <View className="flex-row items-center p-3 border-t border-gray-700">
        <TextInput
          className="flex-1 p-3 bg-gray-800 text-white border border-gray-600 rounded-full mr-2"
          placeholder="Type your message..."
          placeholderTextColor="#a3a3a3"
          value={inputText}
          onChangeText={setInputText}
          onFocus={() => {
            // Stop the typing effect when user starts typing
            setWelcomeMessage('')
            setHelpMessage('')
          }}
        />
        <TouchableOpacity
          onPress={handleUserInput}
          className="p-3 bg-gradient-to-r from-purple-500 to-indigo-600 rounded-full"
        >
          <Ionicons name="send" size={24} color="white" />
        </TouchableOpacity>
      </View>
    </KeyboardAvoidingView>
  )
}

export default ChatScreen
