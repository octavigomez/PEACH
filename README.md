# PEACH-Paleoseismic-EArthquake-CHronologies

This code is designed to compute paleoearthquake chronologies from paleoseismic data, including site chronologies previously modelled with OxCal. 

This work is part of a paper submitted to Earth and Planetary Sciences Letters called: "Deciphering past earthquakes from paleoseismic records – The Paleoseismic EArthquake CHronologies (PEACH) code".

Authors: Octavi Gómez-Novell, Bruno Pace, Francesco Visini, Joanna Faure Walker and Oona Scotti.

Funding: This work has been supported by two consecutive postdoctoral grants awarded to Octavi Gómez-Novell in 2022: “Borsa di studio n. 2341/2021” funded by the INGEO Department (Università degli Studi “G. d’Annunzio” di Chieti e Pescara) and “Margarita Salas grant” funded by the University of Barcelona with EU’s “Next Generation” funds.

Read the user manual in the main folder for details on how to run the code. The folder structure shown in the repository should not be changed.

******

THIRD-PARTY CODES

The principal codes (PEACH.m and PEACH_Oxcal.m) use one external script (allVL1; Bruno Luong, 2022) from the MATLAB Central File Exchange. This code generates all numeric combinations that add to a specified number, and it used to generate all combinations of paleoseismic events that add to the total amount of events in case of overlaps (see user manual).

IMPORTANT. Carefully read the license for the allVL1 code as detailed below:


- allVL1 function: Bruno Luong (2022). All Permutations of integers with sum criteria (https://www.mathworks.com/matlabcentral/fileexchange/17818-all-permutations-of-integers-with-sum-criteria), MATLAB Central File Exchange. Retrieved November 25, 2022.


Copyright (c) 2009, Bruno Luong
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution
* Neither the name of FOGALE nanotech nor the names of its
  contributors may be used to endorse or promote products derived from this
  software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

******

SHARING

The code and derivatives are published under the Creative Commons license CC-BY-NC 4.0. For more info: https://creativecommons.org/licenses/by-nc/4.0/deed.es This means that you are free to copy, share and edit the material as long as you give credit to the authors, publish it under the same license and always for non-comercial purposes. Third party codes (see above) and their licenses must also be disclosed.
