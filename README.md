# MICB475_24W1_Team1

## Tuesday Dec 3 Meeting Agenda/Minutes
1. Look through alpha diversity metrics - should we include Chao1 and Shannon in our manuscript in addition to Faith's?
2. Functional analysis
   <img width="412" alt="Screenshot 2024-12-03 at 10 34 09 AM" src="https://github.com/user-attachments/assets/596eed04-15ac-4f2c-b3cf-470a43ada9c1">


## Tuesday Nov 26 Meeting Agenda/Minutes
1. Go over all analyses and decide on what to include in paper
- ANCOM: <img width="913" alt="Screenshot 2024-11-25 at 2 10 25 PM" src="https://github.com/user-attachments/assets/b0a265da-10a3-46f2-b463-e8f1cc5f7208">
     Conclusions:
     - ART responsive microbiomes more closely resemble healthy microbiomes
     - nonresponsive individuals show further changes in the microbiome after ART/Cotri treatment
 
Minutes/Notes:
* Make diversity plot for the metric that showed how/where the cohorts were needed to be split in R
* Make the weighted unifrac plot in R

- Figure 1
    - Beta diversity metrics (six panels total). Make them only color coded by the response_patient_by_visit so there are no shapes. Run PERMANOVA analysis on all of these
        - Panel A) Unweighted unifrac (emperor) regenerated in R
        - Keep it color-coded by your category to show that it is not dependent on our response category
        - Panel B) Unweighted left cohort
        - Panel C) Unweighted right cohort
            - Top three panels on top of figure
        - Panel DEF) whole, left, and right cohorts in weighted unifrac
- Figure 2
    - Faith’s PD alpha diversity (do kruskal-wallis stat testing)
        - Panel A shows left cohort
        - Panel B shows right cohort
            -  (stacked on top of each other)
            - We see non-responsive have lower diversity, similar to healthy group
- Table 1
    - ANCOM table
        - Make healthy/responsive in rows (ie flip axes)
            - Reference the actual outputs in the GitHub instead of putting it in the supplemental
- Figure 3
    - Core microbiome (make colors the same and normalize)
        - Panel A)
            - Left before (visit 2)
        - Panel B)
            - Left after (visit 3)
        - Panel C)
            - Right before (visit 2)
        - Panel D)
            - Right after (visit 3)
    - Look into the organisms here in the literature to see where/if there is anything of note

- Supplemental figures
    * S1 is Chao1 and Shannon btw left and right cohorts
    * S2 is taxa bar plots, two panels of left and right cohort
        * Squish them into one taxa bar plot by doing averages of the samples per each group
    * S(3?) Indicator taxa table?
    * S(?) could be the heat maps
        * Justifying why we thought this was a valuable avenue to explore

- Keep in mind for the presentation that we should pick important panels for easy interpretation, maybe even do a heat map for the ancom results?
- Put all figures into a single word doc for our last meeting to go thru each one
- Can meet at the end of the week or Monday to discuss our powerpoint and give pointers
  
        
## Tuesday Nov 19 Meeting Agenda/Minutes

1. Go over left/right cohort split analyses
   - Re-run analyses for the separate cohorts?
2. Plan for functional analysis?

NOTES:
- Perform all pairwise comparisons using ANCOM between groups

Option for writing paper
1. Write about both cohorts, accounting for the confounding variable in our discussion, etc
2. Choose one (left or right) cohort later, once all diversity metrics are completed - only add unwanted cohort analyses as supplemental
3. Settle with one cohort for now (left) and write paper based on that, as it appears to have more interesting data. Will need to account for bias in the paper

## Tuesday Nov 5 Meeting Agenda/Minutes
1. Go over analyses
   - DESEQ gives a lot of hits
      - tried ANCOM, not sure about parameters: structural zeros?
      - ALDEx2? Zicoseq?
   - Taxa bar plots do not seem to show much of interest based on our two response columns - is there anything more that can be done?
   - Alpha and Beta Diversity on R
      - currently shannon and chao1, should we do Peliou's?
      - Beta diversity groups to graph: individual or together? all diversity methods?
    
2. Discovered two clustered cohorts when looking at the unweighted unifrac emperor metric; the division is not associated with any of the metadata categories
   - Decided to retry analyses once newly separated cohorts are obtained from teaching team to reduce confounding variables and see how/if our analyses change
   - Most analyses generally did not show much of interest; DESeq showed some differential abundance, but likely influenced by false postives
        - use ANCOM instead for more precision in data
  
3. Decided we will go through with the functional analysis once the analyses are reattempted
   
## Tuesday October 29 Meeting Agenda/Minutes
1. Last week: made phyloseq object
        - includes variations of the response column in case we want to use them
        <img width="767" alt="Screenshot 2024-10-29 at 8 01 46 AM" src="https://github.com/user-attachments/assets/6c2e4a88-4c32-4c25-a5cb-e0260e690f91">
2. This week: run aim 3 analyses
3. Discuss feedback from the proposal document?


## Tuesday October 22 Meeting Agenda/Minutes
1. Splitting response column into 2 columns?
     - assigns response/nonresponse to metadata from both timepoints from same patient -> easier to compare between timepoints downstream for the 2 groups?
     - same for pretreatment viral load
3. correlation analysis for variables other than viral load
     - CD4 count - significant difference
     - cotrimoxazole use, viral load - trend
     - BMI, sex, IL6 level, CRP level, rural/urban - no effect
5. this week: finish Qiime processing and generate phyloseq object

## Tuesday October 15 Meeting Agenda/Minutes
1. Review data wrangling code
2. Ask if we can use 515-806 database for taxonomic analysis
     "Extracted DNA was PCR-amplified with barcoded primers targeting the V4 region of 16S rRNA gene according to the Earth Microbiome Project 16S Illumina Amplicon protocol with the 515F:806R primer constructs"
3. Determine schedule/assignments to complete project proposal -> goal to finish assigned sections by Friday, Oct 18
   - give an update on our progress, if anyone needs help
4. Determine the groupings for diversity/taxonomic analysis
   - responsive/non-responsive (vs healthy?)
   - before and after
5. Ask about title specificity - how many of our goals do we include?

   Minutes:
   - re-label to be more descriptive so that it's ready for publication; HIV low/high --> pre-treatment; negative --> healthy
   - 2 controls: pre-treatment and healthy control
   - can use 515-806 database
   - for data processing --> we do not need to re-label metadata since we will not be using those figures in the manuscript
   - test all groupings
   - We should definitely do functional analysis
   - heatmap --> shows that the data follows a trend that makes sense; higher viral load will end up having more non-responsive individuals
   - Intro and Background --> WHY we are doing the project, what is the significance, what have they found in literature already, how do we build off of that?
   - Hypothesis: need to come up with one; explain what you are going to see and why based on literature
   - Research Objectives: What is the big question, bigger picture
   - difference between approach and overview: approach is detailed (eg. create taxa bar plots) and flowchart is generic (eg. taxonomic analysis)
   - email demux file to Evelyn to check how it was trimmed; do this before diversity analysis
   - Sampling depth: 16,000; not above 39,000 (may have to exclude the experienced non-responsive group)

## Tuesday October 8 Meeting Agenda/Minutes
1. Finalize research question and aims
   - the paper excludes 21 patients with viral load >200 after ART from their analysis
   - microbiome (and immune/inflammatory/general health measures?) signatures associated with poor response to ART
       - 16 ART-naive individuals that have viral load >200 after 24 weeks on ART (ART non-responders): look at alpha/beta diversity and microbial composition before and after 24 week treatment (compare to changes in those that have viral load <200 after 24 weeks on ART.
       - 6 ART-experienced individuals with viral load >200 at baseline: compare with 24 week ART non-responders -> see effect of length of ART treatment?
       - also look at CD4 subtypes (chronic activation, exhaustion) and inflammatory markers in these groups
   - stratify by sex, age -> depends on sample size

Minutes:

Aim 1 – Redefine metadata (and modify manifest file based off of what we have filtered out)
•	First column of sample ID’s
•	Separate each patient between visits 2 and 3
•	Filter out people with only one visit, as well as people on ART to start. 
•	Use visit 2 as a control and visit 3 as an “after” condition. 
•	Make a new column called “response”, where in visit 1 there are possible  “ctrl low”, “ctrl high”, “ctrl healthy (i.e. HIV neg)” categories
•	 In visit 3 have two categories of “responsive” vs “non-responsive”
-	Use histograms to determine highs and lows and responsive vs non-responsive

Aim 2 – correlation analysis btw being responsive or not and start and final viral load (make a matrix w added colors)
- Heat map to visualize correlation (without ctrl healthy): 4 panels

Aim 3 – based on response column, determine if there are general differences in diversity of their gut and taxonomic

-	Taxa bar plot
-	Alpha beta div metrics
-	DeSEQ: for deternining differentially expressed genes.
-	Core microbiome – Venn diagram to compare how many unique versus shared per condition
-	Indicator taxa to see if there are predictive bugs beyond just this dataset
-	Differential abundance analysis

EXTRA: 
Aim 4 – maybe a functional analysis to see what metabolic pathways are activated to see if there are more pathogenic organisms vs commensal organisms, i.e. healthier gut or not


Do research into HIV and ARTs for proposal; look into UJEMI for inspiration for if loops to source code and such. Decide how to define responsiveness, etc., with the histograms or the literature (there was a heat map done on anemia and infancy in UJEMI)
-	We will get a new checklist for the proposal about data wrangling

  
## Tuesday Oct 1 Meeting Agenda/Minutes
1. Brainstorm project 2 research questions:
Potential directions:
HIV dataset #4: effects of cotrimoxazole and ART on microbiota function in males vs females? Take multiple datasets and compare between cohorts/geography?
  - Cohorts: male/female, ART+/-
  - General health measures: BMI, inflammation (IL6, CRP, Cortisol), lymphocyte %, *viral load*, increase/decrease in viral load
    Action Items: look into general health measures -> figure out what counts as healthy or not (recategorize)
    - look into literature, see what variables have been explored
    - Discuss and decide what variables to focus on, plan aims at next meeting
    
IVF dataset: maybe #3 - ion torrent not illumina = big files, not ideal
  
Combining MS/IBD datasets based on similar treatments:
- Ryan dataset: metadata about medications, categorized by class (immunosuppressive...)
- crossover between groups: untreated groups
- not sure how to tie to MS - no tie with inflammation

  
## Monday Sep 23 - prep for team meeting 1
Potential ideas for project 2
1. Comparing gut microbiome in MS patients and ILD patients (halfvarson dataset) -> using existing datasets
   - must have a couple categories in common to compare
   - maybe nutritional correlation with MS
3. Comparing micrbiome of HIV+ mothers and infants -> find 2 datasets
4. Vaginal microbiome relating to IVF and premature birth -> find 2 datasets
5. Use lung microbiome database from ILD patients, compare to vaping/smoking, pre-existing illnesses that makes patients vulnerable to cystic fibrosis --> find our own dataset

### Potential new datasets:
HIV
- https://www.nature.com/articles/s41598-017-06675-1 (Data: https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA354863)
- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4816837/ (Data: https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA307231)
- https://journals.asm.org/doi/10.1128/mbio.00409-23#tabS1 (Data: https://www.ncbi.nlm.nih.gov/bioproject/530161, Metadata: https://journals.asm.org/doi/suppl/10.1128/mbio.00409-23/suppl_file/mbio.00409-23-s0002.tif
- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10837999/  (Data: https://www.ebi.ac.uk/ena/browser/view/PRJEB66206, Metadata: https://zenodo.org/records/8346401)
- Maybe an interesting paper to identify interesting variables: https://journals.asm.org/doi/10.1128/spectrum.00708-21

IVF
- https://www.frontiersin.org/journals/endocrinology/articles/10.3389/fendo.2023.1200002/full#h6 (Data: https://www.ncbi.nlm.nih.gov/Traces/study/?query_key=3&WebEnv=MCID_66f4652bd2376f7b3cb1b872&o=acc_s%3Aa)
- https://www.frontiersin.org/journals/cellular-and-infection-microbiology/articles/10.3389/fcimb.2021.709372/full#h7 (Data: https://www.ncbi.nlm.nih.gov/Traces/study/?query_key=5&WebEnv=MCID_66f4652bd2376f7b3cb1b872&o=acc_s%3Aa)
- https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-021-01184-w#availability-of-data-and-materials (Data: https://www.ncbi.nlm.nih.gov/Traces/study/?query_key=10&WebEnv=MCID_66f3014b8439573273b0cadb&o=acc_s%3Aa)
- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10715154/ (Data: https://www.ncbi.nlm.nih.gov/Traces/study/?query_key=9&WebEnv=MCID_66f4652bd2376f7b3cb1b872&o=acc_s%3Aa)
- https://www.pnas.org/doi/full/10.1073/pnas.1002611107 (Data: https://www.ncbi.nlm.nih.gov/Traces/study/?query_key=11&WebEnv=MCID_66f4652bd2376f7b3cb1b872&o=acc_s%3Aa)
- https://www.frontiersin.org/journals/medicine/articles/10.3389/fmed.2020.00217/full (Data: https://www.ncbi.nlm.nih.gov/Traces/study/?query_key=15&WebEnv=MCID_66f4652bd2376f7b3cb1b872&o=acc_s%3Aa)
