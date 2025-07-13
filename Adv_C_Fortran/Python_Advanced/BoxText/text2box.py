#!/usr/bin/env python3
"""
ASCII Art Text Generator
Converts input text to large ASCII art using box-drawing characters
"""

class ASCIIArtGenerator:
    def __init__(self):
        # Define character patterns using box-drawing characters
        # Each character is 6 rows tall and variable width
        self.patterns = {
            'A': [
                " █████╗ ",
                "██╔══██╗",
                "███████║",
                "██╔══██║",
                "██║  ██║",
                "╚═╝  ╚═╝"
            ],
            'B': [
                "██████╗ ",
                "██╔══██╗",
                "██████╔╝",
                "██╔══██╗",
                "██████╔╝",
                "╚═════╝ "
            ],
            'C': [
                " ██████╗",
                "██╔════╝",
                "██║     ",
                "██║     ",
                "╚██████╗",
                " ╚═════╝"
            ],
            'D': [
                "██████╗ ",
                "██╔══██╗",
                "██║  ██║",
                "██║  ██║",
                "██████╔╝",
                "╚═════╝ "
            ],
            'E': [
                "███████╗",
                "██╔════╝",
                "███████╗",
                "██╔════╝",
                "███████╗",
                "╚══════╝"
            ],
            'F': [
                "███████╗",
                "██╔════╝",
                "███████╗",
                "██╔════╝",
                "██║     ",
                "╚═╝     "
            ],
            'G': [
                " ██████╗ ",
                "██╔════╝ ",
                "██║  ███╗",
                "██║   ██║",
                "╚██████╔╝",
                " ╚═════╝ "
            ],
            'H': [
                "██╗  ██╗",
                "██║  ██║",
                "███████║",
                "██╔══██║",
                "██║  ██║",
                "╚═╝  ╚═╝"
            ],
            'I': [
                "██╗",
                "██║",
                "██║",
                "██║",
                "██║",
                "╚═╝"
            ],
            'J': [
                "     ██╗",
                "     ██║",
                "     ██║",
                "██   ██║",
                "╚██████╔╝",
                " ╚═════╝ "
            ],
            'K': [
                "██╗  ██╗",
                "██║ ██╔╝",
                "█████╔╝ ",
                "██╔═██╗ ",
                "██║  ██╗",
                "╚═╝  ╚═╝"
            ],
            'L': [
                "██╗     ",
                "██║     ",
                "██║     ",
                "██║     ",
                "███████╗",
                "╚══════╝"
            ],
            'M': [
                "███╗   ███╗",
                "████╗ ████║",
                "██╔████╔██║",
                "██║╚██╔╝██║",
                "██║ ╚═╝ ██║",
                "╚═╝     ╚═╝"
            ],
            'N': [
                "███╗   ██╗",
                "████╗  ██║",
                "██╔██╗ ██║",
                "██║╚██╗██║",
                "██║ ╚████║",
                "╚═╝  ╚═══╝"
            ],
            'O': [
                " ██████╗ ",
                "██╔═══██╗",
                "██║   ██║",
                "██║   ██║",
                "╚██████╔╝",
                " ╚═════╝ "
            ],
            'P': [
                "██████╗ ",
                "██╔══██╗",
                "██████╔╝",
                "██╔═══╝ ",
                "██║     ",
                "╚═╝     "
            ],
            'Q': [
                " ██████╗ ",
                "██╔═══██╗",
                "██║   ██║",
                "██║▄▄ ██║",
                "╚██████╔╝",
                " ╚══▀▀═╝ "
            ],
            'R': [
                "██████╗ ",
                "██╔══██╗",
                "██████╔╝",
                "██╔══██╗",
                "██║  ██║",
                "╚═╝  ╚═╝"
            ],
            'S': [
                "███████╗",
                "██╔════╝",
                "███████╗",
                "╚════██║",
                "███████║",
                "╚══════╝"
            ],
            'T': [
                "████████╗",
                "╚══██╔══╝",
                "   ██║   ",
                "   ██║   ",
                "   ██║   ",
                "   ╚═╝   "
            ],
            'U': [
                "██╗   ██╗",
                "██║   ██║",
                "██║   ██║",
                "██║   ██║",
                "╚██████╔╝",
                " ╚═════╝ "
            ],
            'V': [
                "██╗   ██╗",
                "██║   ██║",
                "██║   ██║",
                "╚██╗ ██╔╝",
                " ╚████╔╝ ",
                "  ╚═══╝  "
            ],
            'W': [
                "██╗    ██╗",
                "██║    ██║",
                "██║ █╗ ██║",
                "██║███╗██║",
                "╚███╔███╔╝",
                " ╚══╝╚══╝ "
            ],
            'X': [
                "██╗  ██╗",
                "╚██╗██╔╝",
                " ╚███╔╝ ",
                " ██╔██╗ ",
                "██╔╝ ██╗",
                "╚═╝  ╚═╝"
            ],
            'Y': [
                "██╗   ██╗",
                "╚██╗ ██╔╝",
                " ╚████╔╝ ",
                "  ╚██╔╝  ",
                "   ██║   ",
                "   ╚═╝   "
            ],
            'Z': [
                "███████╗",
                "╚══███╔╝",
                "  ███╔╝ ",
                " ███╔╝  ",
                "███████╗",
                "╚══════╝"
            ],
            '0': [
                " ██████╗ ",
                "██╔═████╗",
                "██║██╔██║",
                "████╔╝██║",
                "╚██████╔╝",
                " ╚═════╝ "
            ],
            '1': [
                " ██╗",
                "███║",
                "╚██║",
                " ██║",
                " ██║",
                " ╚═╝"
            ],
            '2': [
                "██████╗ ",
                "╚════██╗",
                " █████╔╝",
                "██╔═══╝ ",
                "███████╗",
                "╚══════╝"
            ],
            '3': [
                "██████╗ ",
                "╚════██╗",
                " █████╔╝",
                " ╚═══██╗",
                "██████╔╝",
                "╚═════╝ "
            ],
            '4': [
                "██╗  ██╗",
                "██║  ██║",
                "███████║",
                "╚════██║",
                "     ██║",
                "     ╚═╝"
            ],
            '5': [
                "███████╗",
                "██╔════╝",
                "███████╗",
                "╚════██║",
                "███████║",
                "╚══════╝"
            ],
            '6': [
                " ██████╗ ",
                "██╔════╝ ",
                "███████╗ ",
                "██╔═══██╗",
                "╚██████╔╝",
                " ╚═════╝ "
            ],
            '7': [
                "███████╗",
                "╚════██║",
                "    ██╔╝",
                "   ██╔╝ ",
                "   ██║  ",
                "   ╚═╝  "
            ],
            '8': [
                " █████╗ ",
                "██╔══██╗",
                "╚█████╔╝",
                "██╔══██╗",
                "╚█████╔╝",
                " ╚════╝ "
            ],
            '9': [
                " █████╗ ",
                "██╔══██╗",
                "╚██████║",
                " ╚═══██║",
                " █████╔╝",
                " ╚════╝ "
            ],
            ' ': [
                "   ",
                "   ",
                "   ",
                "   ",
                "   ",
                "   "
            ],
            '!': [
                "██╗",
                "██║",
                "██║",
                "╚═╝",
                "██╗",
                "╚═╝"
            ],
            '?': [
                "██████╗ ",
                "╚════██╗",
                "  ▄███╔╝",
                "  ▀▀══╝ ",
                "  ██╗   ",
                "  ╚═╝   "
            ],
            '.': [
                "   ",
                "   ",
                "   ",
                "   ",
                "██╗",
                "╚═╝"
            ],
            ',': [
                "   ",
                "   ",
                "   ",
                "   ",
                "▄██╗",
                "╚═╝"
            ],
            ':': [
                "   ",
                "██╗",
                "╚═╝",
                "██╗",
                "╚═╝",
                "   "
            ],
            ';': [
                "   ",
                "██╗",
                "╚═╝",
                "▄██╗",
                "╚═╝",
                "   "
            ],
            '-': [
                "        ",
                "        ",
                "███████╗",
                "╚══════╝",
                "        ",
                "        "
            ],
            '_': [
                "        ",
                "        ",
                "        ",
                "        ",
                "███████╗",
                "╚══════╝"
            ]
        }
    
    def generate_ascii_art(self, text):
        """Generate ASCII art from input text"""
        text = text.upper()
        
        # Filter out unsupported characters
        supported_chars = []
        for char in text:
            if char in self.patterns:
                supported_chars.append(char)
            elif char == '\n':
                supported_chars.append('\n')
        
        if not supported_chars:
            return "No supported characters found in input text."
        
        # Split text into lines
        lines = ''.join(supported_chars).split('\n')
        result_lines = []
        
        for line in lines:
            if not line:
                # Add empty line
                result_lines.extend([''] * 6)
                continue
            
            # Generate 6 rows for this line
            line_rows = [''] * 6
            
            for char in line:
                if char in self.patterns:
                    pattern = self.patterns[char]
                    for i in range(6):
                        line_rows[i] += pattern[i]
            
            result_lines.extend(line_rows)
        
        return '\n'.join(result_lines)
    
    def generate_with_spacing(self, text, extra_char_spacing=0, extra_line_spacing=0):
        """Generate ASCII art with custom spacing"""
        text = text.upper()
        
        # Split text into lines
        lines = text.split('\n')
        result_lines = []
        
        for line_idx, line in enumerate(lines):
            if not line:
                # Add empty line
                result_lines.extend([''] * (6 + extra_line_spacing))
                continue
            
            # Generate 6 rows for this line
            line_rows = [''] * 6
            
            for char_idx, char in enumerate(line):
                if char in self.patterns:
                    pattern = self.patterns[char]
                    for i in range(6):
                        line_rows[i] += pattern[i]
                        # Add extra spacing between characters (except last char)
                        if char_idx < len(line) - 1:
                            line_rows[i] += ' ' * extra_char_spacing
            
            result_lines.extend(line_rows)
            
            # Add extra line spacing between lines (except last line)
            if line_idx < len(lines) - 1:
                result_lines.extend([''] * extra_line_spacing)
        
        return '\n'.join(result_lines)
    
    def get_supported_characters(self):
        """Return list of supported characters"""
        return sorted(self.patterns.keys())


def main():
    generator = ASCIIArtGenerator()
    
    print("ASCII Art Text Generator")
    print("=" * 50)
    print(f"Supported characters: {', '.join(generator.get_supported_characters())}")
    print()
    
    while True:
        try:
            text = input("Enter text to convert (or 'quit' to exit): ").strip()
            
            if text.lower() == 'quit':
                break
            
            if not text:
                print("Please enter some text.")
                continue
            
            print("\nGenerated ASCII Art:")
            print("-" * 50)
            ascii_art = generator.generate_ascii_art(text)
            print(ascii_art)
            print("-" * 50)
            print()
            
        except KeyboardInterrupt:
            print("\n\nGoodbye!")
            break
        except Exception as e:
            print(f"Error: {e}")


if __name__ == "__main__":
    main()