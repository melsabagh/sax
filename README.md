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

See [docs](https://htmlpreview.github.io/?https://github.com/melsabagh/sax/blob/master/docs/docs.html) for detailed description.

---
### Example: 
```C++
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

