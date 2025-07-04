#!/usr/bin/env python3
# This script creates a FORD project.md from README.md

import re

def remove_title(markdown_text):
    # Remove the project title from the markdown text.

    lines = markdown_text.split('\n')
    lines.remove("# ColSpecF")
    return '\n'.join(lines)

def convert_latex_delimiters(markdown_text):
    # Convert inline LaTeX delimiters from $...$ to \(...\).

    pattern = r'(?<!\$)\$(?!\$)(.*?)(?<!\$)\$(?!\$)'
    return re.sub(pattern, r'\\(\1\\)', markdown_text)

def main(input_file, output_file):
    with open(input_file, 'r', encoding='utf-8') as file:
        markdown_content = file.read()
    
    processed_content = remove_title(markdown_content)
    processed_content = convert_latex_delimiters(processed_content)
    
    with open(output_file, 'w', encoding='utf-8') as file:
        file.write(processed_content)

if __name__ == "__main__":
    main('README.md', 'project.md')
