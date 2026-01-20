# Scientific Considerations & Interpretation Guidance

The Ecosystem Integrity Index represents an ambitious effort to synthesize multiple dimensions of ecosystem health into a single, globally consistent metric. Like any index that attempts to distill complex ecological systems into quantitative scores, the EII involves methodological choices and data constraints that users should understand when interpreting results. This document provides transparent guidance on these considerations.

We present these considerations as important context that shapes how the EII should be applied and interpreted across different use cases. Many of these considerations reflect inherent challenges in global-scale ecological assessment that affect all comparable indices, while others point toward specific areas where ongoing research and data improvements will strengthen future versions.

---

## Understanding the Compositional Pillar

The compositional integrity component draws on the Biodiversity Intactness Index, a metric developed through years of research by the Natural History Museum London and collaborators. The BII estimates how much of an ecosystem's original species community remains relative to an undisturbed reference state; a fundamentally important question for understanding ecosystem health. However, several characteristics of this data source merit consideration when interpreting compositional integrity scores.

### The Challenge of Biodiversity Data

Global biodiversity assessment faces a fundamental data challenge that affects all efforts to map species composition at scale. The statistical models underlying the BII are trained on the PREDICTS database, which represents one of the most comprehensive compilations of biodiversity survey data ever assembled, containing millions of records from thousands of sites worldwide. Nevertheless, the distribution of ecological research effort across the globe is uneven. Biodiversity surveys tend to concentrate in regions with strong research infrastructure, i.e. Western Europe, North America, and well-established tropical research stations, while vast areas of Central Asia, the Sahel, parts of South America, and much of the Middle East have received comparatively less systematic study.

This geographic pattern in underlying data means that the statistical relationships between human pressures and biodiversity responses are calibrated primarily on certain ecosystem types and biogeographic regions. When the BII model predicts compositional integrity in under-sampled regions, it necessarily extrapolates from relationships learned elsewhere. For most applications, this extrapolation is reasonable, namely that the fundamental ecological principle that habitat conversion reduces native species abundance holds globally. However, the magnitude of these relationships may vary across biogeographic provinces in ways that current data cannot fully capture. Users working in regions with sparse biodiversity survey coverage should interpret compositional scores as informative estimates rather than precise measurements, and consider supplementing EII assessments with local biodiversity expertise where available.

### Temporal Considerations

The BII dataset currently used in this implementation extends through 2020-2022. This creates a temporal offset with the functional integrity pillar, which incorporates satellite observations through 2025. For most landscapes, where conditions change gradually, this offset has minimal practical impact. However, in rapidly transforming regions, active deforestation frontiers, areas undergoing major infrastructure development, or sites experiencing rapit restoration interventions, users should recognize that compositional scores may reflect conditions from several years prior. The same is true for the structural connectivity that is estimated from 2022 estimates of human footprint.
Fortunately, the functional integrity pillar, with its near-real-time satellite observations, can serve as an early indicator of ecosystem change, prompting closer investigation in areas where productivity patterns shift unexpectedly.

---

## The Functional Integrity Reference State

Functional integrity compares observed ecosystem productivity to a modeled "natural potential": what productivity we would expect if the landscape were operating under minimal human influence given its climate, soils, and topography. This comparison requires training a predictive model on reference areas that approximate natural conditions, which introduces considerations worth understanding.

### Learning from Earth's Remaining Wild Places

The natural productivity model draws its understanding of ecosystem potential from the planet's most intact landscapes, i.e. areas with minimal human modification, including protected wilderness, remote forests, and other regions where natural processes predominate. This approach is scientifically sound: these areas provide our best available empirical evidence for how ecosystems function absent substantial human alteration. However, the geographic distribution of remaining intact ecosystems is not uniform across all environmental conditions.

In regions with long histories of intensive land use, such as Western Europe, the Indian subcontinent, or eastern China, very few reference areas remain that represent the full range of environmental conditions present in those regions. The productive lowlands and fertile valleys that humans preferentially settled and cultivated are precisely the conditions least represented in remaining natural areas. The model must therefore learn the potential productivity of these environments partly from analogous conditions in other parts of the world, e.g. temperate forests in Eastern North America informing estimates for Western Europe, for instance, or Asian subtropical forests providing context for predictions in southern China.

This extrapolation is guided by well-established ecological principles: productivity responds to temperature, water availability, and soil characteristics in broadly predictable ways across similar environments worldwide. Nevertheless, regional idiosyncrasies in species pools, evolutionary history, and biogeochemical cycling may create local variations that global models cannot fully capture. In heavily transformed regions, functional integrity scores should be interpreted as estimates informed by global ecological relationships rather than locally calibrated benchmarks.

### A Contemporary Reference Point

The EII measures ecosystem condition relative to the best remaining examples of nature that exist today, not relative to some historical pre-industrial baseline. Today's reference landscapes have themselves experienced background changes, elevated atmospheric CO2, nitrogen deposition from industrial emissions, and climate warming, that distinguish them from ecosystems of centuries past. This design choice reflects both practical necessity and a forward-looking philosophy: the EII asks "how does this ecosystem compare to the best we can currently achieve?" rather than "how does it compare to an idealized past we cannot directly observe?"

This framing has important implications for interpretation. An ecosystem that matches contemporary reference conditions receives high functional integrity scores even if historical records suggest that region once supported even higher productivity. The EII is best understood as measuring present-day relative condition rather than absolute deviation from historical baselines.

### Climate Variability and Change

Ecosystem productivity naturally fluctuates with weather patterns: wet years boost plant growth while droughts suppress it. The EII's multi-year averaging of observed productivity (three years) and natural productivity (ten years) dampens the influence of short-term weather variability, focusing the assessment on underlying ecosystem condition rather than transient meteorological events. Additionally, because the natural productivity model incorporates climate variables, its predictions adjust to reflect what should be expected under current climatic conditions.

Nevertheless, climate change creates interpretive nuances worth considering. In regions experiencing persistent climate shifts: prolonged drought in semi-arid zones, altered precipitation patterns, or temperature changes affecting growing seasons, observed productivity may decline for reasons related to regional climate trends rather than local land management. The EII does not currently distinguish between productivity reductions caused by climate change versus those caused by direct human degradation of the land surface. Both represent departures from historical ecosystem function, but they imply different response strategies. Users interpreting declining functional integrity trends in climate-sensitive regions should consider whether regional climate data supports a climate-driven explanation alongside or instead of local degradation hypotheses.

---

## Structural Integrity: Configuration and Its Limits

The structural integrity pillar assesses landscape fragmentation and habitat quality through a core area approach. That means, it is identifying interior habitat patches that are buffered from edge effects and weights them by their degree of human modification. This methodology captures important aspects of landscape ecology while necessarily simplifying others.

### What Core Area Captures

The core area approach reflects well-established ecological principles. Habitat edges experience altered microclimates, elevated predation pressure, invasion by disturbance-adapted species, and other effects that penetrate some distance into habitat patches. Small fragments may consist entirely of edge-affected habitat, while large contiguous areas maintain interior conditions that support species requiring undisturbed environments. By eroding habitat patches to identify core areas, the EII's structural pillar distinguishes between landscapes with equivalent total habitat area but different configurations; recognizing that one large intact patch generally supports greater ecological integrity than many small fragments totaling the same area.

The 300-meter edge depth used in the calculation represents a reasonable central estimate based on ecological literature, where most documented edge effects occur within this distance. However, edge penetration varies considerably across ecosystem types and taxonomic groups. Tropical forests may experience edge effects extending 500 meters or more, while grassland edges may stabilize within 100-150 meters. The globally consistent parameterization chosen for the EII represents a principled compromise that performs reasonably across diverse ecosystems while acknowledging that region-specific edge depths might improve local accuracy.

### The Connectivity Question

Structural integrity as currently implemented measures landscape configuration -- the size, shape, and quality of habitat patches -- but does not directly assess functional connectivity between patches. Two landscapes could receive similar structural integrity scores while differing markedly in their connectivity: two large patches separated by 20 kilometers of inhospitable terrain versus two large patches 500 meters apart across a permeable agricultural matrix. Both landscapes contain equivalent core habitat area, but organisms' ability to move between patches differs substantially.

This distinction matters because connectivity influences population viability, gene flow, and species' capacity to track shifting climates. Measuring functional connectivity rigorously requires species-specific information about dispersal abilities and matrix permeability; data that cannot feasibly be incorporated at global scale with current knowledge. The EII's structural pillar should therefore be understood as measuring habitat configuration and intactness rather than functional connectivity per se. For applications where connectivity is paramount, wildlife corridor planning, climate adaptation assessments, or population viability analyses, users should complement EII structural scores with dedicated connectivity modeling using species-specific parameters.

---

## Relationships Among the Three Pillars

A conceptual strength of the EII framework is its integration of three distinct dimensions of ecosystem health. In an ideal implementation, these pillars would capture independent information, with each contributing unique insight to the overall assessment. In practice, the pillars share some underlying drivers, creating correlations that warrant acknowledgment.

Both structural integrity and compositional integrity are influenced by human pressure data: the Human Modification Index underlies structural calculations directly, while similar pressure variables inform the statistical models predicting biodiversity intactness. This shared causal pathway means these two pillars are not fully independent: landscapes with high human modification tend to score lower on both structural and compositional metrics. The functional pillar, derived from satellite observations of actual vegetation productivity, operates more independently, though it too responds to severe landscape modification.

This correlation does not invalidate the three-pillar framework structure and composition remain conceptually distinct aspects of ecosystem health even if they respond to similar pressures. However, it suggests that the effective dimensionality of the EII may be somewhat less than three fully independent axes. For most applications, this consideration is minor: the framework still captures more ecosystem complexity than single-metric alternatives, and the different pillars respond at different rates to disturbance and recovery, providing complementary temporal signals.

---

## Scale, Resolution, and Boundary Effects

The EII operates at 300-meter resolution globally; a scale that represents careful balancing of competing considerations. Finer resolution would capture more local detail but would dramatically increase computational requirements for global coverage and introduce additional noise. Coarser resolution would smooth over meaningful local variation in ecosystem condition.

At 300 meters, each pixel represents approximately nine hectares; sufficient to characterize landscape units meaningful for many management applications, but too coarse to resolve fine-scale features. Linear infrastructure like roads and power lines, small clearings, and narrow riparian corridors may fall below the resolution threshold. In landscapes with fine-grained heterogeneity, e.g. European agricultural mosaics, peri-urban interfaces, or regions with complex land tenure creating small management units, the 300-meter pixel may blend distinct conditions into intermediate values that don't precisely represent any actual land parcel.

Boundary effects also merit consideration. At transitions between contrasting land covers, e.g. urban edges, forest-agriculture boundaries, coastlines, pixels may span multiple ecosystem types, creating scores that reflect the mixture rather than either component. For applications requiring precise site-level assessment, our proposal of a local modulation layer comes into play [to be documented]. In general, where applicable, aggregating EII scores to property or management unit boundaries often provides more robust characterization than individual pixel values.

---

## Temporal Alignment and Lag Effects

The three EII pillars integrate data from different temporal windows, reflecting the varying update frequencies of their underlying datasets. Functional integrity draws on satellite observations updated continuously, structural integrity uses human modification data updated periodically, and compositional integrity, as chosen for the first round of implementation, relies on biodiversity models with their own update cycles. For most landscapes, where conditions evolve gradually, this temporal spread creates minimal interpretation challenges.

However, in locations undergoing rapid change, the pillars may temporarily reflect different time periods. A recently cleared forest might show collapsed productivity (current satellite data) while structural and compositional scores still reflect pre-clearance conditions from earlier datasets. This temporal asynchrony typically resolves within a few years as all data sources update, but users monitoring dynamic landscapes should recognize that pillar scores may temporarily disagree during transitional periods.

A related consideration involves ecological time lags. Biodiversity responds to habitat change with delays that can span years to decades: the "extinction debt" phenomenon where species persist for some time after their habitat becomes unsuitable, and the "colonization credit" where species gradually accumulate following habitat improvement. The BII's design accounts for expected equilibrium conditions rather than transient states, meaning compositional scores may anticipate changes that have not yet fully manifested in actual species presence. This characteristic is both a strength (providing forward-looking assessment) and a limitation (potentially diverging from current realized conditions during transitions).

---

## The Aggregation Philosophy

The EII combines its three pillars using a "minimum with fuzzy logic" approach, reflecting the ecological principle that ecosystem health is constrained by its most degraded dimension. This aggregation philosophy is deliberately conservative: it prevents an ecosystem with severe degradation in one aspect from receiving high overall scores simply because other aspects remain intact. A productive monoculture plantation, for instance, might score well on functional integrity but would be penalized by low compositional scores, resulting in a modest overall EII.

This design choice embeds a normative assumption that all three dimensions of ecosystem health matter, and that severe deficiency in any dimension represents meaningful degradation regardless of performance in other dimensions. For some applications, this assumption aligns well with user objectives. For others, particularly those focused on specific ecosystem services like carbon sequestration where functional integrity might be paramount, users may find value in examining pillar-level scores directly in addition to the aggregated index.

The fuzzy logic component modulates the minimum-based aggregation, ensuring that ecosystems degraded across multiple dimensions receive lower scores than those degraded in only one. This refinement captures the intuition that cumulative degradation is worse than isolated problems, while maintaining the conservative property that no ecosystem can score higher than its lowest-performing pillar.

---

## Guidance for Different Applications

The EII performs differently across various use cases, and matching the index to appropriate applications maximizes its value while avoiding over-interpretation.

For broad-scale monitoring and reporting, national or continental assessments, corporate portfolio screening, identification of priority regions, the EII provides valuable standardized comparisons across diverse landscapes. At these scales, individual pixel uncertainties average out, and the consistent methodology enables meaningful comparisons that would be difficult to achieve with heterogeneous local assessments.

For regional conservation prioritization and landscape planning, the EII offers useful guidance when interpreted as relative ranking rather than absolute measurement. Comparing EII scores across candidate sites within a region provides defensible prioritization, even if absolute scores carry uncertainty. The pillar decomposition can inform strategy: regions where structural integrity lags might benefit from connectivity restoration, while those with low compositional scores might prioritize habitat protection for remnant biodiversity.

For site-level assessment and project monitoring, the EII is best used as one input among several. The 300-meter resolution, temporal considerations, and inherent uncertainty in global models all suggest that high-stakes site-level decisions warrant local adaptaion, possibly ground-truthing and local expertise alongside EII scores.

For restoration monitoring, users should recognize that the compositional pillar may respond slowly to habitat improvements due to ecological time lags and can, by design, only respond to recovery driven by mappable drivers (primarily threat reduction); but not, for example direct actions, such as, for example species reintroduction or anti-poaching activities. Functional integrity, with its continuous satellite observations, often provides earlier signals of ecosystem recovery. Tracking pillar-level trends separately, rather than relying solely on the aggregated index, can provide more nuanced understanding of restoration trajectories.

---

## Ongoing Development and Future Directions

The considerations outlined in this document inform an active research agenda aimed at strengthening future EII versions. Priorities include expanding the geographic and taxonomic coverage of biodiversity training data, improving temporal alignment across data sources, developing region-specific parameterizations where data support calibration, and implementing uncertainty quantification that communicates confidence levels alongside point estimates.

We are also investigating approaches to better separate climate-driven productivity changes from land-management effects, exploring connectivity metrics that could complement the current configuration-based structural assessment, and working toward more frequent updates of all component layers. User feedback on which improvements would most benefit specific applications is welcome and informs development priorities.

The EII represents our current best effort to provide globally consistent, scientifically grounded assessment of ecosystem integrity. Like all scientific tools, it will evolve as data improve and methods advance. We offer it in a spirit of transparency, acknowledging its limitations while believing it provides meaningful insight for the many users seeking to understand and improve ecosystem conditions worldwide.

---

## References

Gassert, F., Mazzarello, J., & Hyde, S. (2022). Global 100m Projections of Biodiversity Intactness for the years 2017-2020. Impact Observatory Technical White Paper.

Hill, S. L. L., Harrison, M. L. K., Maney, C., et al. (2022). The Ecosystem Integrity Index: a novel measure of terrestrial ecosystem integrity. *bioRxiv*. doi:10.1101/2022.08.21.504707

Hudson, L. N., Newbold, T., Contu, S., et al. (2014). The PREDICTS database: a global database of how local terrestrial biodiversity responds to human impacts. *Ecology and Evolution*, 4(24), 4701-4732.

Kennedy, C. M., et al. (2019). Managing the middle: A shift in conservation priorities based on the global human modification gradient. *Global Change Biology*, 25(3), 811-826.

Miguet, P., et al. (2016). What determines the spatial extent of landscape effects on species? *Landscape Ecology*, 31(6), 1177-1194.

Ries, L., et al. (2004). Ecological responses to habitat edges: mechanisms, models, and variability explained. *Annual Review of Ecology, Evolution, and Systematics*, 35, 491-522.

Scholes, R. J., & Biggs, R. (2005). A biodiversity intactness index. *Nature*, 434(7029), 45-49.

Theobald, D. M., et al. (2020). Earth transformed: detailed mapping of global human modification from 1990 to 2017. *Earth System Science Data*, 12, 1953-1972.
