import { Stack } from 'expo-router'
import '../../global.css'

export default function RootLayout() {
  return (
    <Stack initialRouteName="index">
      <Stack.Screen
        name="index"
        options={{
          title: 'Afya',
          headerStyle: { backgroundColor: 'black' },
          headerTintColor: 'purple',
          headerTitleAlign: 'center',
          headerTitleStyle: { fontWeight: 'bold', fontSize: 28 },
          animation: 'slide_from_bottom',
        }}
      />
    </Stack>
  )
}
