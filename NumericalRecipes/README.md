## Step 1: Compile the Source File into an Object File

Use the gcc compiler to compile your source file into an object file:

```
gcc -c nrutil.c -o nrutil.o
```

## Step 2: Create a Static Library

Create a static library from the object file using the ar command:

```
ar rcs libnrutil.a nrutil.o
```

## Step 3: Move the Linker 

For a system wide usage, move the linker to `/usr/local/lib`

```
sudo mv libnrutil.a /usr/local/lib
```

## Step 4. Using Environment Variables

Set the LD_LIBRARY_PATH environment variable to include the directory where your library is located. This is useful if you don't want to move the library to a standard directory.

```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
```


