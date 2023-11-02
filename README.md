![Turing patterns](https://github.com/ThomasEWoolley/Mathematical_biology/blob/main/Pictures/Spots_to_stripes_2.png)
# Mathematical Biology
These files provide an entire set of notes and problem sheets for a course in Mathematical Biology. The course covers topics such as constitutive laws, steady states, discrete modelling and partial differential equations.

These notes were written by Dr Thomas E. Woolley. Although not based on a specific book, they are most closely aligned to 
Murray, J. D.. Mathematical Biology I. An Introduction. Vol. 17. : Springer New York, 2002. 
Murray, J. D.. Mathematical Biology II: Spatial Models and Biomedical Applications. Vol. 18. : Springer New York, 2003. 
Beyond these they are the culmination of being taught by a diverse range of undergraduate lecturers. As such, input is acknowledged from lecturers Jim D. Murray, Philip K. Maini, Ruth E. Baker and Eamonnn A. Gaffney. Sincere apologies to anyone who has been off this list.
All errors are the fault of Dr Thomas E. Woolley.

# Folders
The contents of the folders should be pretty self explanatory. However, I expand on the titles here.

- Latex

This folder contains the main latex files that produce the notes. The main latex file "Master.tex" is split into chapters (see below). The notes are "gapped", such that some of the text is whited out to be inputted by the lecturer. The coloour of the text can be controlled from the "Master.tex" file in the line \definecolor{COLR}{rgb}{x,x,x}.

Namely, if the x is set to 0 (i.e. \definecolor{COLR}{rgb}{0,0,0}) then the text colour will be black and so the notes will be complete as in "Full.pdf". However, if x is set to 1 (i.e. \definecolor{COLR}{rgb}{1,1,1}) then the text will will be gapped as in "Gapped.pdf".

- Maple

Maple codes that do various algebraic manipulations and Taylor series that are too laborious to be done by hand.

- Matlab

Matlab codes that produce all the figures and animations that are contained within the notes.

- Matlab_students

A subsection of the codes from the Matlab file that are more animation orientated. These codes can be provided to students to play with.

- Pictures

All the figures that appear in the notes and problem sheets.

- Problem_sheets

Five problem sheets, with answers and computational problems covering the main topics of the course.

Like the latex notes the question and answers are provided in a single .tex document. However, the answers can be hidden by changing the line \includecomment{Answ} to \excludecomment{Answ}.

![Simulated diffusion with different boundary conditions](https://github.com/ThomasEWoolley/Mathematical_biology/blob/main/Matlab/Diffusions.gif)
# Course contents
Below I reproduce the sections that appear in the latex notes.

## 1 Introduction
1 Introduction

1.1 References

## 2 Things you have forgotten
2.1 How to model a system

2.2 Law of Mass Action

2.3 Non-dimensionalisation

2.3.1 Examples of non-dimensionalisation through substitution of variables

2.3.2 Examples of non-dimensionalisation through the arrow method

2.4 Stationary states and stability

2.5 Linear stability

2.6 Curve sketching

2.6.1 Specific curves and their properties

2.7 Phase planes

2.7.1 One-dimensional

2.7.2 Two-dimensional

2.8 Check list

## 3 Population modelling
3.1 Continuum modelling

3.1.1 Disease transmission

3.2 Discrete modelling

3.2.1 Cobwebbing

3.2.2 Stability of discrete equations

3.2.3 Oscillatory states in discrete equations

3.3 Check list

## 4 Spatial systems
4.1 Deriving the diffusion equation

4.2 Deriving the taxis equation

4.3 Travelling waves

4.3.1 Fisher’s equation

4.4 Check list

## 5 Pattern Formation
5.1 French flag patterning

5.1.1 Digit specification in the limb bud

5.2 The Turing instability

5.2.1 Diffusion Driven Instability

5.2.2 Stability without diffusion

5.2.3 Instability with diffusion

5.2.4 Corollaries to Turing’s theory

5.2.5 A comment on domain size and spatial dimensions

5.3 Do they exist?

## A Proof of stability criterion for scalar ODE equations
## B Proof of stability criterion for ODE systems
## C Characterising the stability of a two-dimensional ODE system
