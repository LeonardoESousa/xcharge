<h1 align="center">KMC</h1
  
<p align="center">KMC is a platform to perform directional-rate exciton dynamics via Kinectic Monte Carlo (KMC) algorithm for generalized morphologies. </p>

<h4 align="center"> 
	🚧  Under continous development...  🚧
</h4>

[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
[![license](https://img.shields.io/github/license/LeonardoESousa/KMC?style=plastic)]()
[![down](https://img.shields.io/github/downloads/LeonardoESousa/KMC/total?style=plastic)]()
[![maint](https://img.shields.io/maintenance/yes/2021)]()
[![commit](https://img.shields.io/github/last-commit/LeonardoESousa/KMC?style=plastic)]()

Table of Contents
=================
<!--ts-->
   * [Installation](#Installation)
   * [How to Use](#how-to-use)
   * [Basic Use](#basic-use)
   *    * [Make Animations](#ani)
   *    * [Paralellization](#paralell)

<!--te-->

Installation
============

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

How to Use
============

1) Construct your input (look [here](https://github.com/LeonardoESousa/KMC/tree/main/input_examples) for examples)
2) At the same directory of your input, run
```bash
kmc your_input.py
```

Basic Use
============

Make Animations
============
Set, at your input (look [here](https://github.com/LeonardoESousa/KMC/tree/main/input_examples) for examples)

```python
animation_mode = True
save_animation = True
```
Then run the job
```python
kmc your_input.py
```

Paralellization
============
Set the number of rounds as you wish, and assign to n_proc the number of avaliable cores (-1 if you want to use the entire machine). Then, turn off the animation_mode. Your input should look something like
```python
animation_mode = False
rounds = 10000  # repeating the dynamics 10000 times 
n_proc = 30     # using 30 cores to execute the job
```
Then run the job
```python
kmc your_input.py
```

### Autors
---

[![Badge](https://img.shields.io/badge/Email_Leo!-%237159c1?style=for-the-badge&logo=gmail&logoColor=red)](mailto:leonardo.sousa137@gmail.com)

[![Badge](https://img.shields.io/badge/Email_Tiago!-%237159c1?style=for-the-badge&logo=gmail&logoColor=red)](mailto:tiagofis96@gmail.com)
