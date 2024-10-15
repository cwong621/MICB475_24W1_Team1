# MICB475_24W1_Team1

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
