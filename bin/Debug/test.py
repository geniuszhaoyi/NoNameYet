from ctypes import *
dll = cdll.LoadLibrary('dll_test.dll')
dll.test.restype=c_char_p
dll.test.argtype=[c_char_p]
ret = dll.test('{"key":"what is it? ","specie":"E.coli","location":"336..2798","pam":"NGG"}',5)
print ret;
