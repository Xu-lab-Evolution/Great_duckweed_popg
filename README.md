<h1> Population Genomics of the Great Duckweed </h1>

> ⚠ **Important Notice**: The raw WGS sequencing data for 131 <em>Spirodela</em> samples may not be easily found for download on NCBI. Please use the following Run Selector link for access: [NCBI Run Selector](https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA701543)

The repository archives scripts involved in the great duckweed (<strong><em>Spirodela polyrhiza</em></strong>) population genomics project. 
The paper now can be found here: [link](https://www.nature.com/articles/s42003-024-06266-7)



<h2> Prerequisites </h2>
The scripts were written mainly using bash, Perl, Python and R, and calculations were done in an x86_64 GNU/Linux cluster. Please modify the running environment accordingly.

<h2> How to start </h2>
Before running the scripts, please download the raw data and store the "data" folder in the same directory where you put the "scripts" folder, like this:


<ul>
  <li>
    <span class='tree-item'>working_folder/</span>
    <ul>
      <li>
        <span class='tree-item'>scripts/</span>
        <ul>
          <li>
            <span class='tree-item'>Fig_1/</span>
          </li>
          <li>
            <span class='tree-item'>Fig_2/</span>
          </li>
          <li>
            <span class='tree-item'>...</span>
          </li>
        </ul>
      </li>
      <li>
       <span class='tree-item'>data/</span>
       <ul>
         <li>
            <span class='tree-item'>Fig_1/</span>
          </li>
          <li>
            <span class='tree-item'>Fig_2/</span>
          </li>
          <li>
            <span class='tree-item'>...</span>
          </li>
       </ul>
     </li>
    </ul>
  </li>
</ul>
<strong>The tar.gz compressed folder of "data" could be downloaded <a href="https://irods-web.zdv.uni-mainz.de/irods-rest/rest/fileContents/zdv/project/m2_jgu-evoltroph/Duckweed_popg/228_Sp_popg.data_scripts_for_publish/data.tar.gz?ticket=jqq2ixjb6RqNFBg" target="_blank" rel="noopener noreferrer">here</a>.</strong> Or, an alternative link: [10.5281/zenodo.14009270](https://doi.org/10.5281/zenodo.14009270)

<h2> Other files </h2>

Folder <a href="https://github.com/Xu-lab-Evolution/Great_duckweed_popg/tree/main/S.polyrhiza_SpGA2022_annotation" target="_blank" rel="noopener noreferrer"> "S.polyrhiza_SpGA2022_annotation/" </a> contains the updated <em>S. polyrhiza</em> genome annotation "SpGA2022".


Folder <a href="https://github.com/Xu-lab-Evolution/Great_duckweed_popg/tree/main/Supplementary_data" target="_blank" rel="noopener noreferrer"> "Supplementary_data/" </a> contains the supplementary data files involved in the duckweed_popg paper.


