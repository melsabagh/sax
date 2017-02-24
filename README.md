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
#include <string>
#include <iostream>
#include "saxquantizer.hpp"

using namespace std;

int main(int argc, char **argv) {
  SaxQuantizer::Sax sax(6, 1, 5);
  vector<double> input = { 1, 1, 1, 1, 2, 2, 10, 10, 100 };
  vector<int> output;
  sax.quantize(input, &output, false); // true for reduction

  string delim = "";
  cout << "input:";
  for (const auto & x : input) cout << " " << x;

  cout << endl;

  cout << "output:";
  for (const auto & x : output) cout << " " << x;
  cout << endl;
}
```

