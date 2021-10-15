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

<p align="center">
  <kbd>
    <img width="640" style="border-radius: 25px" height="480" src="https://i.postimg.cc/2yVWMfZc/ezgif-com-gif-maker.gif" alt="demo">
  </kbd>
</p>

## Main Features

- [x] Singlet exciton diffusion with Förster transfer rates.
- [x] Triplet exciton diffusion with Dexter and/or Triplet-to-Singlet transfer rates.
- [x] Miller-Abrahams Rate
- [x] Charge dissociation and recombination
- [x] Fully customable morphology (presets: Bulk Heterojunction, Bilayer and more)
- [x] Lattice distortion effects
- [x] Built-in Parallelization
- [x] Singlet Annihilation
- [ ] Include in-site dipole orientation


## Table of Contents

<!--ts-->
   * [Installation](#installation)
   * [How to Use](#how-to-use)
   * [Basic Usage](#basic-usage)
   *    * [Make Animations](#make-animations)
   *    * [Parallelization](#parallelization)
   * [Documentation](#documentation)

<!--te-->

## Installation

(Manual installation)
```bash
git clone https://github.com/LeonardoESousa/KMC
cd KMC/
pip3 install .
```
(Auto installation)
```bash
pip3 install kmc
```

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

Set the number of rounds as you wish and assign to n_proc the number of avaliable cores (-1 if you want to use the entire machine). Then, turn off the animation_mode. Your input should look something like
```python
animation_mode = False
rounds = 10000  # repeating the dynamics 10000 times 
n_proc = 30     # using 30 cores to execute the job
```
Then, run the job
```python
kmc your_input.py
```


## Documentation

<!--ts-->
   * [Input's Basic Structure](#input)
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

