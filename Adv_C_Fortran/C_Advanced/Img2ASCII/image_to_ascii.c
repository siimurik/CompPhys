// Compile and execute with:
//  >gcc -o image_to_ascii image_to_ascii.c -lm
// ./image_to_ascii path_to_your_image_file.jpg

// Define STB_IMAGE_IMPLEMENTATION before including the header to ensure the implementation part
// of the stb_image library is included in this file. This defines the actual functions.
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include <stdio.h>
#include <stdlib.h>

// Function to map grayscale values to ASCII characters
// The function takes a grayscale value (0-255) and returns a corresponding ASCII character.
// The characters are chosen based on their perceived brightness.
char get_ascii_char(int gray_value) {
    const char* ascii_chars = "@%#*+=-:. "; // Characters sorted by brightness
    int num_chars = 10; // Number of characters in the ASCII gradient
    return ascii_chars[(gray_value * (num_chars - 1)) / 255]; // Map grayscale value to ASCII character
}

// Function to convert an image to ASCII art
void image_to_ascii(const char* filename) {
    int width, height, channels;
    // Load the image using stb_image. It returns a pointer to the image data.
    // The function also sets the width, height, and number of channels in the image.
    unsigned char* img = stbi_load(filename, &width, &height, &channels, 0);
    if (img == NULL) {
        // If the image loading fails, print an error message and exit.
        printf("Error in loading the image\n");
        exit(1);
    }

    if (channels < 3) {
        // If the image doesn't have at least 3 channels (RGB), print an error message and free the image memory.
        printf("The image does not have enough color channels\n");
        stbi_image_free(img);
        exit(1);
    }

    // Iterate over the pixels of the image. We downsample vertically by 2 for a better aspect ratio in text.
    for (int y = 0; y < height; y += 2) {
        for (int x = 0; x < width; x++) {
            // Calculate the index in the image data array for the current pixel
            int idx = (y * width + x) * channels;
            int r = img[idx];       // Red component
            int g = img[idx + 1];   // Green component
            int b = img[idx + 2];   // Blue component

            // Convert to grayscale by averaging the RGB values
            int gray_value = (r + g + b) / 3;
            // Get the corresponding ASCII character for the grayscale value
            char ascii_char = get_ascii_char(gray_value);
            printf("%c", ascii_char); // Print the ASCII character
        }
        printf("\n"); // New line after each row
    }

    // Free the image memory allocated by stb_image
    stbi_image_free(img);
}

int main(int argc, char** argv) {
    if (argc < 2) {
        // If the program is called without a filename argument, print usage information and exit.
        printf("Usage: %s <image file>\n", argv[0]);
        return 1;
    }

    // Call the function to convert the provided image file to ASCII art
    image_to_ascii(argv[1]);
    return 0;
}
