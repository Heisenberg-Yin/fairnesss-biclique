import os
from collections import defaultdict


with open('output_2.txt') as file_obj:
  lines = file_obj.readlines()

line_dict=defaultdict(int)  
for line in lines:
  line_dict[line]=1
with open('output.txt') as file_obj:
  lines = file_obj.readlines()
for line in lines:
  if(line_dict[line]==0):
    print(line)