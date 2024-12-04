# Phytochemical diversity and seasonality drive a critical transition in the structure of a plant-herbivore network

This README file was generated on 14.08.2024 by Leandro Giacobelli Cosmo

GENERAL INFORMATION

1. Title of Dataset: 

Cosmo_et_al_critical_transition_network

2. Author information:

Leandro G. Cosmo: Programa de Pós-Graduação em Ecologia, Departamento de Ecologia, Instituto de Biociências, Universidade de São Paulo, São Paulo, Brazil / Department of Evolutionary Biology and Environmental Studies, University of Zurich, Winterthurerstrasse 190, Zurich CH-8057, Switzerland

Kate P. Maia: Departamento de Ecologia, Instituto de Biociências, Universidade de São Paulo, São Paulo, Brazil.

Paulo R. Guimaraes Jr.: Departamento de Ecologia, Instituto de Biociências, Universidade de São Paulo, São Paulo, Brazil.

Martin Pareja: Departamento de Biologia Animal, Instituto de Biologia, Universidade Estadual de Campinas, Campinas, Brazil.

Corresponding author: Leandro G. Cosmo, E-Mail: leandro.giacobellicosmo@uzh.ch

2. Information about funding sources that supported the collection of the data:

This study was financed in part by the Coordenação de Aperfeiçoamento de Pessoal de Nível Superior – Brasil (CAPES) – Finance Code 001. LGC was funded by a São Paulo Research Foundation PhD scholarship (FAPESP; grants # 2019/22146-3 and #2022/07939-0), and is currently funded by a SNSF postdoctoral fellowship (Swiss National Science Foundation, grant #197201). KPM is funded by a FAPESP postdoctoral felloswhip (grant #2019/21732-6). PRG is funded by CNPq (307134/2017-2), FAPESP (São Paulo Research Foundation; grant: # 2018/14809-0), and the Royal Society, London (CHL/R1/180156).

DATA & FILE OVERVIEW

1. File List: 

Functions:

criteria.R: Function to compute the connectivity parameter as defined in the main text (equation 1).

katz_function: Function to compute the number of direct and indirect pathways in the network, as defined in the main text (equation 3).

summarySE: Function to compute summary statistics.

Datasets:

network_dry.txt/network_rainy.txt: Matrix of interactions between plant individuals (rows) and herbivore species (columns) at the dry and rainy seasons of the study site, respectively. Each entry of the matrix corresponds to the abundance of a given herbivore species found feeding on a given plant individual.

plant_dataset.csv: Dataset containing information of each plant individual used on the study for the following variables. (1) Identification of the plant individual (plant); (2) season of data collection (season); (3) percentage of leaf damaged (herbivory); (4) structural phytochemical diversity, as defined in the main text (structural_pd); and (5) compositional phytochemical diversity, as defined in the main text (compositional_pd).

Scripts:

script_main_analyses.R: Script used to perform the main analyses in the main text and to plot figures 2-4.

script_sensitivity_analysis_bootstrap.R: Script used to perform the sensitivity analysis as defined in the Supplementary Information.
