<h1 align="center">Xcharge: A Kinetic Monte Carlo Model for Exciton and Charge Dynamics</h1
  
<p align="center">Xcharge is a Kinectic Monte Carlo (KMC) algorithm that can be used to study exciton (singlet or triplet) and charge (electron or hole) dynamics in various morphologies. </p>

<h4 align="center"> 
	🚧  Under continuous development...  🚧
</h4>

[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
[![license](https://img.shields.io/github/license/LeonardoESousa/KMC?style=plastic)]()
[![down](https://img.shields.io/github/downloads/LeonardoESousa/KMC/total?style=plastic)]()
[![maint](https://img.shields.io/maintenance/yes/2021)]()
[![commit](https://img.shields.io/github/last-commit/LeonardoESousa/KMC?style=plastic)]()

<!--
<p align="center">
  <kbd>
    <img width="640" style="border-radius: 25px" height="480" src="https://i.postimg.cc/2yVWMfZc/ezgif-com-gif-maker.gif" alt="demo">
  </kbd>
</p>
-->

<!--
<p align="center">
  <kbd>
    <img width="441" style="border-radius: 25px" height="340" src="https://i.postimg.cc/Zqnz3gnS/final.gif" alt="demo">
  </kbd>
  &nbsp;&nbsp;&nbsp;&nbsp;
  <kbd>
    <img width="407" style="border-radius: 25px" height="357" src="https://i.postimg.cc/s2hZVnzV/name-cut.gif" alt="demo">
  </kbd>
</p>
-->

<p align="center">
  <kbd>
    <!--<img width="399.3" style="border-radius: 25px" height="308.55" src="https://i.postimg.cc/Zqnz3gnS/final.gif" alt="demo">-->
    <!--<img width="399.3" style="border-radius: 25px" height="308.55" src="https://i.postimg.cc/7P7KWFJM/kmc.gif" alt="demo">-->
    <img width="399.3" style="border-radius: 25px" height="308.55" src="https://i.postimg.cc/C1yPf4t0/kmc-opt.gif" alt="demo">
  </kbd>
  &nbsp;&nbsp;&nbsp;&nbsp;
  <kbd>
    <!--<img width="335.5" style="border-radius: 25px" height="293.7" src="https://i.postimg.cc/s2hZVnzV/name-cut.gif" alt="demo">-->
    <img width="335.5" style="border-radius: 25px" height="293.7" src="https://i.postimg.cc/zfNHbbjM/latt.gif" alt="demo">
  </kbd>
</p>



## Main Features

- [x] Easy to use!
- [x] Singlet exciton diffusion with Förster transfer rates;
- [x] Triplet exciton diffusion with Dexter and/or Triplet-to-Singlet transfer rates;
- [x] Miller-Abrahams Rate;
- [x] Charge dissociation and recombination;
- [x] Fully customable morphology (presets: Bulk Heterojunction, Bilayer and more);
- [x] Lattice distortion effects;
- [x] Built-in Parallelization;
- [x] Singlet Annihilation;
- [x] Native Dashboard for Analysis; 
- [ ] Include in-site dipole orientation.


## Table of Contents

<!--ts-->
   * [Installation](#installation)
   * [How to Use](#how-to-use)
   * [Basic Usage](#basic-usage)
   *    * [Make Animations](#make-animations)
   *    * [Parallelization](#parallelization)
   *    * [Analysis Dashboard](#analysis-dashboards)
   * [Documentation](#documentation)

<!--te-->

## Installation

(Manual installation)
```bash
git clone https://github.com/LeonardoESousa/KMC
cd KMC/
pip3 install .
```
<!--
(Auto installation)
```bash
pip3 install kmc
```
-->

If you encounter any error, we recommend you to visit our page of [troubleshooting solutions](https://github.com/LeonardoESousa/xcharge/wiki/Troubleshoot-Install)
## How to Use

1) Construct your input (look [here](https://github.com/LeonardoESousa/KMC/tree/main/input_examples) for examples)
2) At the same directory, run
```bash
kmc your_input.py
```

## Basic Usage

### Make Animations

Set, at your input (look [here](https://github.com/LeonardoESousa/KMC/tree/main/input_examples) for examples),

```python
animation_mode = True
save_animation = True
```
Then run the job
```python
kmc your_input.py
```

### Parallelization

Set the number of rounds as you wish and assign to n_proc the number of avaliable cores. Then, turn off the animation_mode. Your input should look something like
```python
animation_mode = False
rounds = 10000  # repeating the dynamics 10000 times 
n_proc = 30     # using 30 cores to execute the job
```
Then, run the job
```python
kmc your_input.py
```

### Analysis Dashboard

<p align="center">
  <kbd>
    <img width="640" style="border-radius: 25px" height="480" src="https://i.postimg.cc/fRvGx1sK/dash.gif" alt="demo">
  </kbd>
</p>

At any terminal, enter
```bash
dash
```
Click in <em> Upload </em>, choose you file and press <em>Read File </em>.

## Documentation
[Main page](https://github.com/LeonardoESousa/xcharge/wiki)
<!--ts-->
   * [Input's Basic Structure](https://github.com/LeonardoESousa/xcharge/wiki/Input's-Basic-Structure)
   * [Morphology Design](#morph)
   * [Energy Setting](#energy)
   * [Particle Creation](#particles)
   * [Rates](#rates)
<!--te-->

## Authors
---
<table>
  <tr>	  
    <td align="center"><a href="https://github.com/LeonardoESousa"><img style="border-radius: 50%;" src="https://avatars.githubusercontent.com/u/49243510?v=4" width="100px;" alt=""/><br /><sub><b>Dr. Leonardo Evaristo de Sousa</b></sub></a><br />Postdoctoral researcher at the Techinical University of Denmark (DTU)<br /><sub><b><a href="mailto:leonardo.sousa137@gmail.com">EMAIL</a><br /></td>
    <td align="center"><a href="https://github.com/TSA-Cassiano"><img style="border-radius: 50%;" src="https://avatars.githubusercontent.com/u/26448170?s=400&u=b0820613fd46515f0cfe1806f7251a414d4c249b&v=4" width="100px;" alt=""/><br /><sub><b>MSc Tiago de Sousa Araújo Cassiano</b></sub></a><br />Phd student at the University of Brasília (UnB)<br /><sub><b><a href="mailto:tiagofis96@gmail.com">EMAIL</a><br /></td>
  </tr>
</table>

