SVMAsm
======

Spam filter based on SVM, written in C++ and assembly language(NASM).


## Compiling

For building application, the NASM assembler is required. Moreover, the algorithms are implemented in 32-bit assembly, so you need 32-bit GNU libraries on 64-bit systems. For example, on Debian-based distros the mandatory package is:
* g++-multilib

After that, you can just type
```
cd src/ && make
```

## Testing

Directory `test` contains unit tests, written with the gtest framework, and test data for learning.
