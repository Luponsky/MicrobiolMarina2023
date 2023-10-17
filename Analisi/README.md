### SEA4BLUE Project

<p align="center">
  <img width="460" height="300" src="https://github.com/Luponsky/MicrobiolMarina2023/blob/main/img/Sea4Blue_page-0001.jpg">
</p>

Sail4Blue samples were collected from Florida to the Azores during May and June 2022 along the transect illustrated in the figure above. The samples were collected using the eDNA Citizen Scientist Sampler with the Self-preserving eDNA Filter. Each sample corresponds to 8 to 10 liters of filtered water. The filters were stored at room temperature (RT), and DNA was extracted using the DNeasy PowerWater kit by Qiagen.


<p align="center">
  <img width="502" height="206" src="https://github.com/Luponsky/MicrobiolMarina2023/blob/main/img/gel.png">
</p>
The 16S rRNA V4 region was amplified starting with 5 ng of each sample 


## The data we will use has already been 'pre-processed' with DADA2 and the taxonomy was assigned

<details>
<summary>What are the 'pre-processing' steps needed to clean the raw reads?</summary>
<br>
Inspect read quality profiles<br>
Remove Adaptors<br>
Quality filter and trimming<br>
Sample Inference<br>
Merge paired reads<br>
Remove chimeras<br>
</details>


<details>
<summary>What are the outputs of DADA2?</summary>
<br>
ASVs table & Rapresentative Sequences
</details>


