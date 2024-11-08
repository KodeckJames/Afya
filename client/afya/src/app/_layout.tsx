import { Stack } from 'expo-router'
import '../../global.css'

export default function RootLayout() {
  return (
    <Stack initialRouteName="index">
      <Stack.Screen
        name="index"
        options={{
          title: 'Afya',
          headerStyle: { backgroundColor: '#1F2937' },
          headerTintColor: '#B794F4',
          headerTitleAlign: 'center',
          headerTitleStyle: { fontWeight: 'bold', fontSize: 28 },
          animation: 'slide_from_bottom',
        }}
      />
    </Stack>
  )
}
