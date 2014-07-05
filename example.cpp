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


