**Reproducible Figures Are More Valuable Than Fancy Figures**

**Reproducibility**  
Reproducibility, a fundamental characteristic of legitimate scientific methods (Popper, 1959), refers to the ability of independent investigators to achieve consistent (or comparable) outcomes upon replicating an experiment or analysis. It represents an indispensable element of the scientific methodological framework. In recent years, both researchers and academic journals have actively sought to ensure reproducibility by proactively making research data, materials, and analytical code accessible (Stodden, Seiler, & Ma, 2018). This practice provides scientists with evidence supporting the credibility and authenticity of research findings, thereby mitigating the influence of bias or uncontrolled chance (Gannot _et al_., 2016). Consequently, reproducible results are amenable to independent verification, critical evaluation, and further development, thereby reinforcing the foundation of scientific knowledge (Peng, 2011).

**Irreproducibility and the Reproducibility Crisis**  
According to Landis et al. (2012) and Shamoo and Resnik (2015), irreproducibility, on the other hand, might indicate a problem with any of the steps involved in the research, including the experimental design, variability of biological materials (such as cells, tissues, or animal or human subjects), data quality or integrity, statistical analysis, study description, or incomplete methods. As numerous studies have revealed that a sizable percentage of published research findings cannot be repeated, questions regarding irreproducibility have grown in recent years (Ioannidis, 2005; Baker, 2015; Collins & Tabak, 2014). The "reproducibility crisis" in science is a term frequently used to describe this increasing awareness.

**Scientific Misconduct and Retractions**  
Data fabrication or falsification may contribute to some of the irreproducibility in scientific research (Collins & Tabak, 2014; Kornfeld & Titus, 2016; Shamoo, 2013, 2016). Over two-thirds of retractions are the result of negligence or alleged misconduct (Fang, Steen, & Casadevall, 2012). Papers that have been withdrawn because of misbehavior or other issues are methodically tracked by the website Retraction Watch (Retraction Watch, 2016).

**Reproducibility in Computational and Data-Driven Research**  
According to Peng (2011) and Sandve, Nekrutenko, Taylor, and Hovig (2013), reproducibility in computation-based research means that results derived from data must be supported by raw data and code that would allow the study to be repeated from start using established techniques. To put it another way, repeatability plays a major role in determining how transparent the data is (Stodden et al., 2016). If a scientific study's analysis codes, procedures, data, and findings are accessible for anybody to verify and comprehend, it can be considered reliable and repeatable (Peng, 2011).

**Challenges of Irreproducible Computational Outputs**  
Concerns over the reproducibility of research findings have emerged as experimental outcomes in a variety of scientific domains have increasingly failed to be repeated (Baker, 2015). Instead of reliable data and figures, results are frequently produced as files and fancy images (Sandve et al., 2013). These findings are typically shown or kept in notebooks. Because notebooks are software, they can cause problems when interacting with other software (Rule et al., 2019). Due to their ease of use and need for programming language proficiency, notebooks are primarily utilized by researchers.

**Reproducibility Issues in R-Based Analyses**  
According to a study conducted on 2000 datasets written in R, only 26% of the data were repeatable without hitting mistakes at the time of final publication (Ana et al, 2022). For reproducibility to be satisfied, codes must be executed again (Peng, 2011). In R, developing reproducible artifacts entails gathering R package versions as well as system dependencies used in the code, ensuring that all data sources are included, and closely examining all potential sources of errors and non-conformance (Gentleman & Temple Lang, 2007).

**Reproducibility as an Ethical Issue**  
Reproducibility is an ethical as well as a scientific concern (Shamoo & Resnik, 2015). Scientists may suspect data fabrication or falsification when they are unable to replicate a study result (Kornfeld & Titus, 2016). Reproducibility problems have resulted in accusations of data falsification or fabrication in several well-known situations. In 1986, post-doctoral researcher Margot O'Toole accused her supervisor, assistant professor of pathology at Tufts University Thereza Imanishi-Kari, of fabricating and falsifying data in a study funded by the National Institutes of Health (NIH) on the use of foreign genes to stimulate the production of antibodies in mice, which was published in the journal _Cell_ (Kevles, 1998). When O'Toole failed to replicate a crucial experiment carried out by Imanishi-Kari and discovered differences between the data documented in Imanishi-Kari's lab notebooks and the data presented in the article, she started to have doubts about the study. Concerns regarding irreproducibility have led to inquiries into scientific misconduct in several well-known examples, highlighting the significance of transparent procedures and repeatable analyses in preserving public confidence in research (Shamoo & Resnik, 2015).

**Figures and Their Role in Scientific Communication**  
The purpose of scientific figures is to communicate significant patterns and insights that might not be apparent in the raw data (Torres _et al_., 2025). However, an overemphasis on visual appeal often obscures the primary goal of creating the figures. According to Sandve et al. (2013), figures that have been manually altered or "polished" using external tools may appear professional, but they lack transparency and reproducibility without the underlying source disclosed.

**Script-Driven Figures and Reproducible Visualization**  
R-generated script-driven figures provide a completely distinct method. They can be precisely recreated, examined during peer evaluation, and systematically altered when new data or parameters are introduced because they are created directly from data using documented code (Wickham, 2016). In contrast to manual figure modification, parameterized plots enable researchers to update analyses effectively while preserving consistency among figures (Sandve et al., 2013).

**Reproducibility Over Aesthetics**  
Reproducibility should never be sacrificed for aesthetic refinement, even when visual clarity is crucial. Not only is a figure impossible to replicate, but it is also scientifically fragile (Peng, 2011). Reproducible figures, on the other hand, promote verification, encourage collaboration, and preserve the integrity of scientific research (Stodden et al., 2016).

**References**

Ana Trisovic, Matthew K Lau, Thomas Pasquier, and Mercè Crosas. 2022. A large-scale study on research code quality and execution. _Scientific Data_ 9, 1 (2022), 60.

Baker, M. (2016). 1,500 scientists lift the lid on reproducibility.

Collins, F. S., & Tabak, L. A. (2014). Policy: NIH plans to enhance reproducibility. _Nature_, _505_(7485), 612-613.

Fang, F. C., Steen, R. G., & Casadevall, A. (2012). Misconduct accounts for the majority of retracted scientific publications. _Proceedings of the National Academy of Sciences_, _109_(42), 17028-17033.

Gannot, G., Cutting, M. A., Fischer, D. J., & Hsu, L. J. (2016). Reproducibility and transparency in biomedical sciences. _Oral diseases_, _23_(7), 813.

Gentleman, R., & Temple Lang, D. (2007). Statistical analyses and reproducible research. _Journal of Computational and Graphical Statistics_, _16_(1), 1-23.

Hwang, S. Y., Yon, D. K., Lee, S. W., Kim, M. S., Kim, J. Y., Smith, L., ... & Ioannidis, J. P. (2023). Causes for retraction in the biomedical literature: a systematic review of studies of retraction notices. _Journal of Korean medical science_, _38_(41).

Ioannidis, J. P. (2005). Why most published research findings are false. _PLoS medicine_, _2_(8), e124.

Kevles, D. J. (2016). _The Baltimore case: A trial of politics, science, and character_. WW Norton & Company.

Kornfeld, D. S. (2019). Research misconduct, NSF v NIH: Its nature and prevalence and the impact of their respective methods of investigation and adjudication. _Accountability in Research_, _26_(6), 369-378.

Landis, S. C., Amara, S. G., Asadullah, K., Austin, C. P., Blumenstein, R., Bradley, E. W., ... & Silberberg, S. D. (2012). A call for transparent reporting to optimize the predictive value of preclinical research. _Nature_, _490_(7419), 187-191.

Peng, R. D. (2011). Reproducible research in computational science. _Science_, _334_(6060), 1226-1227.

Rule, A., Birmingham, A., Zuniga, C., Altintas, I., Huang, S. C., Knight, R., ... & Rose, P. W. (2019). Ten simple rules for writing and sharing computational analyses in Jupyter Notebooks. _PLoS computational biology_, _15_(7), e1007007.

Sandve, G. K., Nekrutenko, A., Taylor, J., & Hovig, E. (2013). Ten simple rules for reproducible computational research. _PLoS computational biology_, _9_(10), e1003285.

Shamoo, A. E., & Resnik, D. B. (2009). _Responsible conduct of research_. Oxford University Press.

Stodden, V., Seiler, J., & Ma, Z. (2018). An empirical analysis of journal policy effectiveness for computational reproducibility. _Proceedings of the National Academy of Sciences_, _115_(11), 2584-2589.

Torres, H., Ozturk, E., Fang, Z., Zhang, N., Cai, S., Sarkar, N., & Coskun, A. F. (2025). What is a "Good" figure: Scoring of biomedical data visualization. _Plos one_, _20_(11), e0336917.

Wickham, H. (2016). Data analysis. In ggplot2: elegant graphics for data analysis (pp. 189-201). Cham: Springer international publishing.
