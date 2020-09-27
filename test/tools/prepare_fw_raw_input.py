#!/usr/bin/env python

from __future__ import print_function
import sys
import re
import six

def main():
  if len(sys.argv) < 2:
    print('Usage: %s text-file' % sys.argv[0])
    return

  # Read lines from file
  filename = sys.argv[1]

  with open(filename) as f:
    lines = f.readlines()

  # Parse lines
  keep = False

  cached_lines = {}

  re_header0 = re.compile('==== Endcap (\d+) Sector (\d+) Hits ====')
  re_header1 = re.compile('==== Endcap (\d+) Sector (\d+) Tracks ====')
  re_content0 = re.compile('\w+ ' * 11 + '\w+')
  re_content1 = re.compile('12345')

  for line in lines:
    m = re_header0.match(line)
    if m:
      # Start keeping lines
      keep = True
      endcap, sector = m.group(1, 2)
      es = (int(endcap) - 1) * 6 + (int(sector) - 1)
      cached_lines[es] = []

    elif re_header1.match(line):
      # Finished keeping lines
      keep = False

      # Find the best sector (highest num of lines)
      best_es = -1
      best_es_len = 0
      for k, v in six.iteritems(cached_lines):
        if best_es_len < len(v):
          best_es_len = len(v)
          best_es = k
      if best_es != -1:
        print('\n'.join(cached_lines[best_es]))

      # Clear lines before the next event
      cached_lines.clear()

    if keep:
      if (re_content0.match(line) or re_content1.match(line)) and not line.startswith('bx'):
        cached_lines[es].append(line.strip())


# ______________________________________________________________________________
if __name__ == '__main__':

  main()
