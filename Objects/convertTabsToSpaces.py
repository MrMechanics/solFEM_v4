# -*- coding: utf-8 -*-
"""
Created on Sat Oct 14 15:48:10 2023

@author: chatGPT
"""

def convert_tabs_to_spaces(filename):
    with open(filename, 'r') as file:
        content = file.read()

    updated_content = content.replace('\t', '    ')

    with open(filename, 'w') as file:
        file.write(updated_content)

if __name__ == '__main__':
    filename = input("Enter the python file name: ")
    convert_tabs_to_spaces(filename)
    print(f"Converted tabs to spaces in {filename}")