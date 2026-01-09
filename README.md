# Filter Diagonalization - Work in Progress

## Overview

This project is a **work in progress** aimed at implementing **filter diagonalization** for quantum chemistry applications. **Filter diagonalization** is a numerical method used to compute eigenvalues for very large matrices, and while it is not meant to replace existing diagonalization methods, it has been very helpful in specific situations where other methods do not perform well.

The method has been particularly useful when coupled with a parallel implementation of my program, which employs a **low storage** approach. Instead of storing the entire matrix, it calculates the needed product at each iteration. This makes it possible to handle **large systems** that would otherwise be impractical to process with conventional methods.

## Key Features

- **Filter Diagonalization**: This project includes a basic implementation of the filter diagonalization method, which is particularly useful when paired with efficient storage methods like low storage techniques, allowing for better performance on large systems.
  
- **Quantum Chemistry Application**: The working version of the program uses **core-ADC (Core-Approximate Diagonalization)**, a quantum chemistry method designed to compute **core-ionization energies** for molecules or clusters. This method has been used successfully at the **University of Perugia** and **Heidelberg University**.

- **C Library**: The filter diagonalization library is written in **C** by a beginner (myself!) back in 2002/2003. While the code is **not elegant** and might not follow best practices, the **date** of the files is significant to me, as they represent my early steps in programming. I was deeply touched when I saw the output-header with my name printed after so many years!

- **Testing and Updates**: The `try/` directory contains tests that use the library. While the code still doesn't work perfectly, it serves as a basis for further improvement. My programming experience has improved since I originally wrote the code, and the new version of the tests reflects this.

## Current Status

- The project is still **in progress**, and there are known issues with the current implementation.
- The filter diagonalization method is not a replacement for other, more established diagonalization methods. In fact, those methods tend to work better when they are applicable.
- The filter diagonalization approach was primarily useful to me due to its compatibility with a **low-storage, parallelized version** of my program, which allowed me to handle very large systems that wouldn't have been feasible otherwise.

## License

This project is distributed under the **[MIT License](LICENSE)**. See the LICENSE file for more details.

## Acknowledgments

- **Prof. Francesco Tarantelli**: This project is dedicated to my PhD supervisor, who first introduced me to Unix computer systems and software development. His enthusiasm for theoretical chemistry deeply influenced my work.
- **Prof. Robin Santra**: I would also like to thank Prof. Santra for suggesting the use of filter diagonalization to extract selected eigenvalues for very large matrices.
- **Prof. Lorenz Cederbaum**: Special thanks to Prof. Cederbaum's research group in Heidelberg, where I had the privilege to work with many brilliant scientists.

## To-Do

- **Improve stability**: The current version of the code still has stability issues and will need further work.
- **Refactor code**: Refactor the existing C code to follow better programming practices.
- **Improve documentation**: Add more documentation and examples to make the code more accessible to others.

## Contact

For any questions or contributions, feel free to reach out to me via **[c.villani@freenet.de]**.

