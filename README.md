This is the overall repository that replicates all analyses from
_Urai AE, Braun A, Donner TH. (2017) Pupil-linked arousal is driven by decision uncertainty and alters serial choice bias. Nature Communications, 8: 14637.
http://www.nature.com/articles/ncomms14637_

### This repository contains ###

- Task_2IFC_RandomDots, containing all the code needed to run the behavioural and eye-tracking task. This code runs under PsychToolbox-3 (http://psychtoolbox.org). The scripts s1a and s1b show instructions and examples to the participant (in Dutch), and s3_2IFC_Main.m runs the main experiment.
- Analysis, containing all Matlab and R files to replicate the analyses and figures in the paper.
- serial-dependencies, a modification of https://bitbucket.org/mackelab/serial_decision/ that includes a modulatory interaction term in the history model.

###What to do to replicate the full analyses###

1. Get the data and code
  1. download and unzip the data files from FigShare _http://dx.doi.org/10.6084/m9.figshare.4300043_, and put everything in a folder called Data. Note that you can skip the individual P??.zip files if you're not interested in redoing the pupil preprocessing and motion filtering. In that case, getting just the GrandAverage.zip and CSV.zip files will do.
  2. get the GitHub code repository and change the name of the folder to Code.
  5. put both Code and Data in a folder called pupilUncertainty, you'll need the name of this path.

2. Rerun the analyses
  * open a0_reproduceAnalyses.m, specify the path to code + data and follow instructions
  * if you're missing functions, they can almost certainty be found in my Tools folder https://github.com/anne-urai/Tools

If you have any questions, open an issue or get in touch @AnneEUrai / anne.urai@gmail.com

###LICENSE###

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
If you use the Software for your own research, cite the paper.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
