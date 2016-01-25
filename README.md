# pupilUncertainty #

This is the overall repository that replicates all analyses from
_Urai AE, Braun A, Donner THD. Decision Uncertainty Drives Pupil-linked Arousal Systems and Modulates Sequential Choice Bias._

### This repository contains ###
- Task_2IFC_RandomDots, containing all the code neede to run the behavioural and eye-tracking task. This code runs under PsychToolbox-3 (http://psychtoolbox.org). The scripts s1a and s1b show instructions and examples to the participant (in Dutch), and s3_2IFC_Main.m runs the main experiment.
The data that were collected are in the Zenodo repository XXX.
- Analysis, containg all Matlab files to replicate the analyses. Going through the files in alphabetical order will get you from data to final figures. The Data folder contains a few sample files to try out the code.
- serial-dependencies, a modification of https://bitbucket.org/mackelab/serial_decision/ that includes a modulatory interaction term in the history model. 


###What to do to replicate our full analyses###

1. Get the data and code
  1. download and unzip all the files from FigShare
  2. put all the folders P01-P27 inside the Data folder to recreate the folder structure we want to work with
  3. put the MotionEnergy folder inside the Data folder
  4. get the repository with the code from GitHub and change the name of the folder to Code
  5. put both Code and Data in a folder called pupilUncertainty, you'll need the name of this path.

2. Rerun the analyses
  * open a0_reproduceAnalyses.m and follow instructions
  * if you're missing functions, they can almost certainty be found in my Tools folder https://github.com/anne-urai/Tools

###LICENSE###

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
If you use the Software for your own research, cite the paper.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
