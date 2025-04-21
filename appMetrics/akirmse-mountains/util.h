/*
 * MIT License
 * 
 * Copyright (c) 2017 Andrew Kirmse
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef _UTIL_H__
#define _UTIL_H__

#include <string>
#include <vector>
#include <map>
#include <unordered_set>
#include <sstream>


/*
 * Split the given string by the given delimiter, putting the
 * result in elems and returning it.
 */
inline std::vector<std::string> &split(const std::string &s, char delim,
                                std::vector<std::string> &elems) {
    elems.clear();
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
      elems.push_back(item);
    }
    return elems;
  }


// Remove the first instance of any (key, value) mapping from the given multimap
template <typename K, typename V>
void removeFromMultimap(std::multimap<K, V> *mmap, K key, V value) {
  const auto range = mmap->equal_range(key);
  auto it = range.first;
  for (; it != range.second; ++it) {
    if (it->second == value) {
      mmap->erase(it);
      break;
    }
  }
}

// Remove all elements of vec whose indices appear in the given set.
template <typename T>
void removeVectorElementsByIndices(std::vector<T> *vec, const std::unordered_set<int> &indices) {
  int to = 0;
  for (int from = 0; from < (int) vec->size(); ++from) {
    if (indices.find(from) == indices.end()) {
      if (from != to) {
        (*vec)[to] = (*vec)[from];
      }
      to += 1;
    }
  }

  vec->erase(vec->begin() + to, vec->end());  // adjust size
}

#endif  // ifndef _UTIL_H__
