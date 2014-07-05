#SaxQuantizer 

A time/space optimized implementation of Symbolic Aggregate approXimation (SAX), with numerosity reduction and scaling. 
This is based on ``A Symbolic Representation of Time Series, with Implications for Streaming Algorithms,'' by Jessica Len et al.

---
### Dependencies

* C++ compiler with C++11 support
* Boost Math

---
### Usage

Simply include `saxquantizer.hpp` in your project, and use the provided functions from the SaxQuantizer namespace. 

---
### Example: 
```
#include <vector>
#include <iostream>
#include "saxquantizer.hpp"

int main(int argc, char **argv) {

  SaxQuantizer::Sax sax(6, 1, 5);
  vector<double> input = { 1, 1, 1, 1, 2, 2, 10, 10, 100 };
  vector<int> output;
  sax.quantize(input, &output, false); // true for reduction

  std::string delim = "";
  std::cout << "input: ";
  for (auto x : input) {
    std::cout << delim << x;
    delim = ", ";
  }
  std::cout << std::endl;

  delim = "";
  std::cout << "output: ";
  for (auto x : output) {
    std::cout << delim << x;
    delim = ", ";
  }
  std::cout << std::endl;

}
```

---
### Public Member Functions

* Sax (size_t window_size, size_t string_size, size_t alphabet_size)
* template<typename Container > void train (const Container &samples)
* template<typename Container > size_t quantize (const Container &seq, vector< int > *qseq, bool reduce=true)
* size_t order () const
* double ratio () const

See `docs` for detailed description.
