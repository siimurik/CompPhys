## Step 1: Compile the Source File into an Object File

Use the `gcc` compiler to compile your source file into an object file:

```bash
gcc -c nrutil.c -o nrutil.o
```

## Step 2: Create a Static Library

Create a static library from the object file using the `ar` command:

```bash
ar rcs libnrutil.a nrutil.o
```

## Step 3: Move the Library and Header Files

For system-wide usage, move the static library and the header file to their appropriate directories:

```bash
sudo mv libnrutil.a /usr/local/lib
sudo mv nrutil.h /usr/local/include
```

## Step 4: Using Environment Variables (Optional)

If you prefer not to move the library and header files to standard directories, you can set environment variables to include their locations. However, this is not necessary if you place them in the standard directories.

### For Library Files:

Set the `LIBRARY_PATH` environment variable to include the library directory:

```bash
export LIBRARY_PATH=$LIBRARY_PATH:/usr/local/lib
```

### For Header Files:

Set the `C_INCLUDE_PATH` environment variable to include the header file directory:

```bash
export C_INCLUDE_PATH=$C_INCLUDE_PATH:/usr/local/include
```

## Step 5: Compile Your Program

Once the files are in the correct locations, you can compile your program without additional flags:

```bash
gcc polint.c romberg.c -lnrutil -lm -o rom
```

