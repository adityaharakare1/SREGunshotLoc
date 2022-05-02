# import socket
# _ipAddr = '169.254.70.92'
# _cnnctPort = 9999
#
# sckt = socket.socket(socket.AF_INET,socket.SOCK_STREAM)
# # sckt.bind((_ipAddr,_cnnctPort))
# # sckt.bind((socket.gethostname(),_cnnctPort))
# sckt.bind(('localhost',_cnnctPort))
# # sckt.listen(2)
# # conn, addr = sckt.accept()
# # while(1):
# #     msg = sckt.recv(1024)
# #     print (msg.decode('utf-8'))

# Import socket module
import socket

# Create a socket object
s = socket.socket()

# Define the port on which you want to connect
port = 12345

# connect to the server on local computer
s.connect(('127.0.0.1', port))

# receive data from the server and decoding to get the string.
print(s.recv(1024).decode())
# close the connection
s.close()