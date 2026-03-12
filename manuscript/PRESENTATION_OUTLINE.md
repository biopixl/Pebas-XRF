# Pebas-XRF Presentation Outline
## Paleoenvironmental Reconstruction from XRF Core Scanning

### Recommended Slide Organization (18 slides)

**Structure Overview:**
| Part | Topic | Slides | Purpose |
|------|-------|--------|---------|
| 1 | Introduction | 1-3 | Context and research questions |
| 2 | Instrumentation & Approach | 4-6 | How XRF works, element/proxy selection |
| 3 | Data Quality | 7 | Core coverage and QC |
| 4 | Results: Site Contrast | 8-9 | TAM vs SC fundamental differences |
| 5 | Validation | 10-11 | Paleontological confirmation |
| 6 | Key Intervals | 12-15 | Specific events and environmental modes |
| 7 | Synthesis | 16-18 | Depositional model and conclusions |

---

## PART 1: INTRODUCTION (3-4 slides)

### Slide 1: Title & Context
- **Title**: "Redox and Carbonate Dynamics in the Miocene Pebas Mega-Wetland: XRF Evidence from Tamshiyacu and Santa Corina"
- Key visual: Map of Western Amazonia showing TAM and SC localities
- Time period: Middle Miocene (~12 Ma), Mollusc Zones MZ7-MZ8

**TALKING POINTS:**
> "We're looking at one of the most remarkable freshwater ecosystems in Earth history - a mega-wetland larger than the modern Mediterranean Sea that persisted for nearly 15 million years."

> "Our cores come from two localities about 30 km apart - close geographically, but as we'll see, representing fundamentally different environments within this vast system."

> "The Tamshiyacu site - TAM - corresponds to the well-studied IQ26 locality where Wesselingh and colleagues documented exceptional mollusc diversity."

### Slide 2: The Pebas System
- >1 million km² freshwater mega-wetland
- Endemic mollusc radiations: 158 species, 72% endemic
- Key references: Wesselingh et al. 2002; Hoorn et al. 2010
- Visual: Pebas paleogeographic reconstruction

**TALKING POINTS:**
> "The Pebas system wasn't a simple lake - it was a complex mosaic of shallow lakes, swamps, embayments, and river channels, all connected in a mega-wetland that covered most of western Amazonia."

> "This system produced one of the most spectacular freshwater mollusc radiations known - 158 species, 72% found nowhere else. For comparison, the African Great Lakes, famous for their cichlid radiations, have far fewer endemic molluscs."

> "The key question is: what environmental conditions allowed this remarkable biodiversity? The mollusc ecology tells us about habitat preferences, but geochemistry can reveal the underlying water chemistry."

**BACKGROUND - Why Pebas matters:**
- Andean uplift created a foreland basin that trapped water
- Marine connections via Venezuela brought occasional brackish influence
- System terminated ~8 Ma when drainage reversed to form the modern Amazon

### Slide 3: Research Questions
1. Can XRF proxies detect the dysoxic conditions documented by mollusc ecology?
2. Do geochemical signals match independent isotopic evidence from bivalve shells?
3. What drove environmental variability in the Pebas system?

**TALKING POINTS:**
> "Wesselingh showed that Pachydon obliquus - the most abundant bivalve - thrives in dysoxic, organic-rich muds. Can we see this signal in the sediment chemistry?"

> "Kaandorp analyzed stable isotopes in the bivalve shells themselves. Do the water chemistry signals from shells match what we see in bulk sediment?"

> "The ultimate question: was this system static, or did it experience dramatic environmental swings? And if so, what drove them?"

---

## PART 2: INSTRUMENTATION & ANALYTICAL APPROACH (3-4 slides)

### Slide 4: XRF Core Scanning - How It Works
**Visual: ITRAX schematic or photo of core scanner**

| Component | Specification | Significance |
|-----------|---------------|--------------|
| X-ray source | Mo tube, 30 kV, 45 mA | Excites K-lines for elements Al through Sr |
| Detector | Si-drift detector (SDD) | Energy-dispersive, resolves individual elements |
| Step size | 2 mm | Each measurement = ~3-5 years at Pebas sed rates |
| Integration | 10 seconds/point | Adequate counts for trace elements |
| Output | Counts per second (cps) | Semi-quantitative, relative abundances |

**TALKING POINTS:**
> "XRF core scanning is non-destructive - the core stays intact. An X-ray beam hits the sediment surface, exciting atoms that fluoresce at characteristic energies. Each element has a fingerprint."

> "We use a molybdenum X-ray tube because Mo excitation efficiently produces K-shell fluorescence for the elements we care about: manganese, iron, calcium, titanium. Lighter elements like aluminum are harder to detect - they're absorbed by air and water."

> "The 2mm step size is a compromise. Finer resolution means more measurements but lower counts per point. At 2mm, we get excellent signal-to-noise while still capturing millimeter-scale environmental changes."

> "Important caveat: XRF gives us counts, not concentrations. We report ratios to cancel out matrix effects and water content variations. The absolute numbers are semi-quantitative; the patterns are robust."

**BACKGROUND - XRF physics:**
- X-ray fluorescence: incident photons eject inner-shell electrons; outer electrons drop down, emitting characteristic X-rays
- Energy-dispersive detection: SDD measures photon energies, software identifies elements by peak position
- Matrix effects: heavier elements absorb X-rays from lighter elements; water attenuates signal
- Why ratios work: matrix effects and geometry cancel when dividing element by element

### Slide 5: Element Selection - Why These Elements?

**Elements measured and their geochemical behavior:**

| Element | Source | Behavior | Use in This Study |
|---------|--------|----------|-------------------|
| **Mn** | Authigenic | Redox-sensitive: precipitates when O₂ present, diffuses when absent | Primary redox indicator |
| **Fe** | Detrital + authigenic | Less redox-sensitive than Mn; forms sulfides under anoxia | Supporting redox indicator |
| **Ca** | Biogenic (shells) | Redox-insensitive; tracks carbonate production/preservation | Productivity/shell proxy |
| **Ti** | Detrital only | Chemically inert; locked in rutile, ilmenite | Normalizing element |
| **Al** | Detrital (clays) | Inert but poorly detected by XRF in wet sediment | Alternative normalizer (not used) |
| **K, Rb** | Detrital (clays, feldspars) | Inert, track clay mineralogy | Considered but not primary |

**TALKING POINTS:**
> "Element selection isn't arbitrary. We need elements that respond to the environmental variable we're interested in - oxygen - and elements that don't respond, to use as normalizers."

> "Manganese is the ideal redox tracer because its solubility changes dramatically at the oxic-anoxic boundary. Mn⁴⁺ in MnO₂ is insoluble - it stays put. Mn²⁺ is soluble - it diffuses into the water and is lost from the sediment."

> "Titanium is nearly indestructible. It sits in heavy minerals like rutile and ilmenite that survive weathering, transport, and diagenesis. Ti doesn't care about oxygen, pH, or salinity - it just tracks how much detrital sediment arrived."

> "We considered aluminum as a normalizer - it's the classic choice in marine geochemistry. But Al detection by XRF is poor in wet sediments because the low-energy Al X-rays are absorbed before reaching the detector. Ti is a better choice for core scanning."

**BACKGROUND - Why not other elements?**
- **Sr**: Tracks carbonate but also substitutes for Ca in shells; not independent of Ca
- **S**: Would indicate sulfide formation under anoxia, but S detection by XRF is poor
- **Zr**: Excellent detrital proxy (in zircon) but grain-size sensitive; coarse zircon creates spikes
- **Si**: Dominated by quartz; detection issues similar to Al

### Slide 6: Proxy Rationale - Building the Interpretive Framework

| Proxy | Numerator Source | Denominator Source | Independence | Interpretation |
|-------|------------------|-------------------|--------------|----------------|
| **Mn/Ti** | Authigenic Mn oxides | Detrital Ti minerals | ✓ Primary | Oxic = high, Dysoxic = low |
| **Ca/Ti** | Biogenic carbonate | Detrital Ti minerals | ✓ Independent | Shells/productivity = high |
| **Fe/Mn** | Mixed Fe sources | Authigenic Mn | ✗ Shares Mn | Supporting only; high = reducing |
| **Fe/Ti** | Mixed Fe sources | Detrital Ti minerals | ✓ Independent | Detrital + authigenic Fe input |

**TALKING POINTS:**
> "The key insight from compositional data analysis: ratios that share an element aren't statistically independent. Fe/Mn and Mn/Ti both contain Mn - if Mn drops, both ratios change even if Fe and Ti are constant."

> "This is called the closure effect, described by Aitchison in 1986. It creates spurious correlations that can fool you into thinking two processes are linked when they're actually measuring the same thing."

> "Our solution: Mn/Ti as the primary redox proxy, Ca/Ti as the independent productivity check. Ca and Mn have completely different source mechanisms - biogenic versus authigenic. When both move together, it's real."

> "We report Fe/Mn for comparison with older literature that uses this ratio, but we interpret it cautiously. High Fe/Mn supports a reducing interpretation, but it's not independent evidence."

**BACKGROUND - Compositional data constraints (Aitchison 1986):**
- Closure: proportions must sum to 100%, creating negative correlations
- Spurious correlation: shared denominators or numerators create mathematical relationships
- Solution: use log-ratios (centered log-ratio, additive log-ratio) for formal statistics
- Pragmatic approach: select ratios with independent source terms; interpret shared-element ratios cautiously

---

## PART 3: DATA QUALITY & CORE COVERAGE (1 slide)

### Slide 7: Core Sections & Stratigraphic Coverage
- **TAM**: 876 measurements across 11 sections (IQ26 site, MZ7-MZ8)
- **SC**: 1,062 measurements across 16 sections
- Quality control: excluded foam-filled zones, QC-failed measurements
- Visual: Section locations on stratigraphic column with mollusc zone boundaries

**TALKING POINTS:**
> "We collected nearly 2,000 XRF measurements at 2mm resolution - that's submillimeter-scale environmental reconstruction, far finer than traditional bulk sampling."

> "The foam-filled zones you see in some cores are where friable sediment collapsed during storage. We excluded these because XRF can't penetrate air gaps and foam - it would measure the foam chemistry, not the sediment."

> "Quality control matters: the ITRAX scanner flags measurements where X-ray counts are too low or the detector saturates. About 3% of our data was excluded by these automatic quality filters."

> "The stratigraphic coverage spans mollusc zones MZ7 and MZ8 - this is the peak of Pebas biodiversity, when the endemic radiation reached its maximum diversity. We're asking: what were the environmental conditions during this remarkable evolutionary event?"

**BACKGROUND - Data quality considerations:**
- Foam-filled zones: friable organic-rich clay collapsed during storage; excluded from analysis
- QC-failed measurements: low counts (<1000 cps total) or detector saturation; automatic exclusion
- Edge effects: first/last 5mm of each section may show X-ray scatter artifacts
- Replicate scans: selected sections scanned twice; reproducibility CV <5% for major elements

---

## PART 4: RESULTS - SITE-LEVEL CONTRAST (2-3 slides)

**NARRATIVE TRANSITION:**
> "Now that we understand the instrumentation and proxy framework, let's examine what the data reveal. The first and most striking result is the fundamental geochemical difference between our two sites - a difference that tells us about oxygen availability in this ancient ecosystem."

### Slide 8: TAM vs SC - The Fundamental Difference
**Key statistics to display:**

| Parameter | TAM | SC | Ratio | Interpretation |
|-----------|-----|-----|-------|----------------|
| Mn/Ti median | 0.18 | 0.37 | 2.1× | SC more oxic |
| Ca/Ti median | 5.04 | 2.94 | 1.7× | TAM more carbonate |
| Ti median (cps) | 9,510 | 13,587 | 1.4× | SC more detrital |
| Mn median (cps) | 1,531 | 4,804 | 3.1× | SC more Mn precipitation |

**Visual**: Side-by-side density plots of Mn/Ti distribution

**TALKING POINTS:**
> "This single slide captures the fundamental environmental difference between our two sites. TAM has half the Mn/Ti ratio of SC - that's not a subtle difference, it's a completely different oxygen regime."

> "Look at the absolute Mn counts: SC has three times more manganese. That Mn had to come from somewhere - it precipitated from oxygenated water. TAM never had enough dissolved oxygen to trap that manganese."

> "But notice TAM has higher carbonate. More shells accumulated there despite - or perhaps because of - the low oxygen. Pachydon thrived in these conditions."

> "The higher Ti at SC tells us this site received more detrital input - it was closer to sediment sources, probably fluvial channels bringing sand and silt from the Andes."

**BACKGROUND - What drives these differences:**
- TAM (IQ26) = distal lacustrine setting, thermally stratified, bottom waters persistently dysoxic
- SC = proximal fluvio-lacustrine, closer to river inputs, seasonally ventilated
- Same mega-wetland, completely different local environments separated by just 30 km

### Slide 9: Site Contrast - Core Image Comparison
**Use Shiny app to generate:**
- TAM-1-2-3B-A (most dynamic TAM section)
- SC-5-6-7ABC-B (best SC section)
- Display side-by-side with Mn/Ti, Ca/Ti, core image

**TALKING POINTS:**
> "The visual difference is striking - TAM shows persistently low Mn/Ti in that blue-gray clay, while SC oscillates between oxic peaks and reducing troughs."

> "Look at the color of the sediment itself. The dark gray to black intervals at TAM correspond to low Mn/Ti - that's organic-rich, reducing mud. The lighter tan intervals at SC correspond to oxic peaks where Mn precipitated."

> "This is core scanning at its best: you can literally see the chemistry in the sediment color, and the XRF confirms it quantitatively."

**BACKGROUND - Reading core lithology:**
- Dark gray/black = high organic carbon, reducing conditions (sulfides may be present)
- Light gray/tan = oxidized, lower organic content
- Shell-rich horizons appear as white bands with elevated Ca
- Color changes often correlate with Mn/Ti transitions at mm-scale

---

## PART 5: VALIDATION WITH PALEONTOLOGICAL EVIDENCE (2-3 slides)

**NARRATIVE TRANSITION:**
> "The site contrast is clear in the geochemistry. But is it real? Does it match what we know from 20 years of paleontological research on these same sediments? This is where independent validation becomes crucial."

### Slide 10: Pachydon obliquus - The Dysoxia Indicator
**Key quote** (Wesselingh 2006a):
> "*Pachydon obliquus* dominates the faunas of dysoxic, organic-rich, clayish intervals"

**XRF validation**:
- TAM sections with Mn/Ti < 0.15 (strongly dysoxic):
  - TAM-1-2-3B-A: 243-444 mm, mean Mn/Ti = 0.078, Fe/Mn = 436
  - TAM-1-2-3B-B: 898-1207 mm, mean Mn/Ti = 0.091, Fe/Mn = 277

**Visual**: Core photo with Pachydon shells + XRF profile showing low Mn/Ti

**TALKING POINTS:**
> "This is the critical validation. Wesselingh studied these bivalves for decades and showed that Pachydon obliquus - the dominant species - is a dysoxia specialist. It lives where other molluscs can't survive."

> "How does Pachydon tolerate low oxygen? It has physiological adaptations: a large siphon to sample water from just above the sediment-water interface, and possibly anaerobic metabolic pathways. It's the bivalve equivalent of a deep-sea specialist."

> "Our XRF data confirms this independently. Where Pachydon dominates, Mn/Ti drops below 0.15. The Fe/Mn ratio exceeds 400 - that's extreme reducing conditions where almost all Mn has diffused away."

> "This is what validation looks like: an ecological signal from shells matches a geochemical signal from mud, using completely independent methods."

**BACKGROUND - Pachydon ecology:**
- Pachydontinae are endemic to Pebas - found nowhere else in the fossil record
- Thick shells suggest calcification in supersaturated waters (high alkalinity)
- Dominance in dysoxic facies suggests competitive exclusion of other taxa
- Sister taxa (Anodontites) occupy normoxic riverine habitats today

### Slide 11: Bivalve Shell Isotopes Confirm XRF Patterns
**From Kaandorp et al. (2005):**

| Bivalve Group | Isotope Pattern | XRF Equivalent |
|---------------|-----------------|----------------|
| Pachydon (lacustrine) | Low amplitude, irregular δ¹⁸O | Stable, low Mn/Ti at TAM |
| Diplodon (fluvial) | High amplitude, seasonal cycles | Variable Mn/Ti at SC |

**Interpretation**: Stratified, poorly-mixed waters at TAM (low Mn/Ti) vs. seasonally-ventilated waters at SC (variable Mn/Ti)

**TALKING POINTS:**
> "Kaandorp measured oxygen isotopes directly in the bivalve shells - this is the water chemistry recorded by the animals themselves during growth. It's completely independent of our sediment geochemistry."

> "Pachydon shells show low-amplitude, irregular δ¹⁸O - no clear seasonal cycle. Why? Because stratified lakes buffer seasonal changes. The bottom water where Pachydon lives doesn't mix with surface water, so it stays chemically stable."

> "Diplodon shows the opposite: high-amplitude seasonal cycles in δ¹⁸O. This is the signature of fluvial systems that respond to seasonal rainfall and dry/wet cycles. The water chemistry changes dramatically through the year."

> "Our XRF tells the same story: stable, low Mn/Ti at TAM matches the dampened isotope signal; variable Mn/Ti at SC matches the dynamic isotope signal. Two independent proxies, same environmental interpretation."

**BACKGROUND - Isotope interpretation:**
- δ¹⁸O in shells tracks water δ¹⁸O and temperature during growth
- Seasonal δ¹⁸O cycles require: (1) seasonal temperature change, OR (2) seasonal water source change
- Stratified lakes dampen both signals in bottom waters
- Sr isotopes in same shells distinguish Andean (high ⁸⁷Sr/⁸⁶Sr) vs. cratonic (low ⁸⁷Sr/⁸⁶Sr) water sources

---

## PART 6: KEY INTERVALS - THE STORY IN DETAIL (4-5 slides)

**NARRATIVE TRANSITION:**
> "Now let's dive into specific intervals that tell the most compelling stories. These are the moments when the Pebas system did something dramatic - and the XRF captured it."

### Slide 12: Dramatic Oxygenation Events
**Highlight these intervals (use Shiny app for visuals):**

| Section | Depth | Mn/Ti Change | Interpretation |
|---------|-------|--------------|----------------|
| **SC-1ABC-2-3C-A-RUN1** | 367 mm | +833% | Most extreme oxygenation pulse |
| **TAM-5AB-6-7-B-RUN2** | 648 mm | +482% | Abrupt ventilation event |
| **SC-5-6-7ABC-B** | 520 mm | Mn/Ti = 4.67 | Highest sustained oxygenation |

**TALKING POINTS:**
> "An 833% increase in Mn/Ti over just a few millimeters - that's one of the most dramatic oxygenation events I've seen in any lacustrine record. Something fundamentally changed in this system."

> "What could cause such rapid oxygenation? Three possibilities: first, a major flood pulse from the Andes bringing cold, oxygen-rich water. Second, lake-level drop exposing bottom sediments to the atmosphere. Third, wind-driven mixing breaking stratification."

> "The sharpness of these transitions is key. Gradual oxygen changes would produce gradual Mn/Ti gradients as the redox front moves through sediment. These are knife-edge transitions - the environment flipped like a switch."

> "At SC, we see these events repeatedly. At TAM, they're rare and less extreme. This reinforces the site contrast: SC was frequently perturbed, TAM remained stubbornly stratified."

**BACKGROUND - Ventilation mechanisms:**
- Andean rivers carry 8-10 mg/L dissolved O₂; stagnant lake bottoms may have <0.5 mg/L
- Cold underflows: dense, oxygen-rich river water can plunge beneath warm lake water
- Wind mixing depth scales with lake fetch and wind speed; large lakes resist mixing
- Modern analogs: Lake Tanganyika shows similar episodic ventilation of anoxic deep waters

### Slide 13: Shell-Rich Horizons (Coquinas)
**Highlight these intervals:**

| Section | Depth | Ca/Ti | Mn/Ti | Environment |
|---------|-------|-------|-------|-------------|
| **SC-5-6-7ABC-B** | 704 mm | 48.4 | 0.50 | Shell lag, moderate O₂ |
| **TAM-1-2-3B-B** | 1204-1213 mm | 34-42 | 0.13-0.20 | Shell bed, dysoxic |
| **SC8-A** | 360-363 mm | 31-32 | 0.23-0.50 | Mixed conditions |

**TALKING POINTS:**
> "Ca/Ti values above 30 are extraordinary - that's 30 times more calcium than titanium. You're essentially looking at carbonate with a little mud mixed in, not mud with some shells."

> "The key insight: high Ca/Ti under low Mn/Ti. How can you have abundant shells in dysoxic water where most animals can't live? Three explanations."

> "First: transported shells. Storm events or density currents wash shells from oxic shallows into the anoxic basin. The shells accumulate, but Mn stays low because the bottom water remains reducing."

> "Second: Pachydon can live there. These dysoxia-tolerant bivalves grow in situ, producing shells even in low-oxygen mud. The shells are autochthonous, not transported."

> "Third: authigenic carbonate. High alkalinity in stratified waters can precipitate calcite directly from the water column. This would explain carbonate without shells."

**BACKGROUND - Shell taphonomy in dysoxic settings:**
- Dysoxia actually favors shell preservation - no bioturbation, reduced dissolution
- Transported shells show abrasion, fragmentation, mixed size classes
- In situ Pachydon beds show articulated valves, life-position orientations
- Coquinas (shell beds) may represent: storm lags, hardgrounds, or winnowed residues

### Slide 14: Concordant Oxic-Carbonate Events (Strongest Signal)
**These intervals show INDEPENDENT proxy agreement:**

| Section | Depth | Mn/Ti (z-score) | Ca/Ti (z-score) | Interpretation |
|---------|-------|-----------------|-----------------|----------------|
| **SC-3AB-4ABCD-RUN2-F** | 1258 mm | +5.8σ | +2.0σ | Textbook oxic-productive event |
| **SC-5-6-7ABC-C** | 848 mm | +2.3σ | +4.4σ | Strong concordance |
| **TAM-5AB-6-7-B-RUN2** | 264 mm | +4.3σ | +2.1σ | Best TAM concordant event |
| **TAM-1-2-3B-A** | 519-522 mm | +3.2-3.4σ | +2.0-2.3σ | Sustained favorable conditions |

**TALKING POINTS:**
> "This is where compositional data analysis pays off. Because Ca and Mn aren't chemically linked, when both ratios move together, it's real signal, not mathematical artifact."

> "A +5.8 sigma anomaly means this Mn/Ti value is nearly 6 standard deviations above the mean - an event that should occur less than once in a million measurements by chance. This is statistically extraordinary."

> "What does concordant Mn/Ti and Ca/Ti mean biologically? Oxygenation favored mollusc growth and shell production. The system switched to a favorable state: well-oxygenated water, productive, shell-producing conditions."

> "These are your highlight intervals for any presentation. They're statistically robust, mechanistically interpretable, and they tell a clear story: the Pebas system experienced discrete favorable episodes in an otherwise challenging environment."

**BACKGROUND - Why concordance matters:**
- Mn and Ca have different source terms: Mn = authigenic precipitation, Ca = biogenic production
- Concordant anomalies require a common environmental driver affecting both processes
- Most parsimonious interpretation: oxygenation events that (1) precipitated Mn and (2) allowed mollusc colonization
- Z-scores normalize each proxy to its own variance, enabling comparison across different scales

### Slide 15: Environmental Modes Summary
**Three distinct modes identified:**

1. **Oxic-Carbonate Events** (Mn/Ti↑, Ca/Ti↑)
   - Enhanced ventilation + productivity
   - Best examples: SC-3AB-4ABCD-RUN2-F @ 1258mm; TAM-5AB-6-7-B-RUN2 @ 264mm

2. **Dysoxic Baseline** (Mn/Ti low, Ca/Ti variable)
   - Persistent reducing conditions, Pachydon habitat
   - Best examples: TAM-1-2-3B-A @ 243-444mm; TAM-1-2-3B-B @ 898-1207mm

3. **Shell Transport/Authigenic Events** (Ca/Ti↑, Mn/Ti low)
   - Carbonate accumulation under reducing conditions
   - Best examples: TAM-1-2-3B-B @ 1204mm (Ca/Ti=42, Fe/Mn=167)

**TALKING POINTS:**
> "These three modes aren't arbitrary categories - they emerge from the data and make ecological sense. The Pebas system oscillated between these states."

> "Mode 1 - the oxic-carbonate state - is rare but spectacular. These are the golden intervals where everything aligned: oxygen, productivity, shell growth. They're the exception, not the rule."

> "Mode 2 - dysoxic baseline - is the dominant state, especially at TAM. This is Pachydon's world: low oxygen, high alkalinity, organic-rich muds. Challenging for most life, but Pachydon thrived."

> "Mode 3 is the puzzle. High carbonate under persistent reducing conditions requires either transport mechanisms or in-situ Pachydon production. The taphonomy of the shells would distinguish these."

> "Notice the modes aren't site-specific. Both sites experience all three modes, but at different frequencies. TAM spends more time in Mode 2; SC oscillates more between Modes 1 and 2."

**BACKGROUND - Environmental mode interpretation:**
- Mode transitions may be driven by: climate (precession cycles), tectonics (Andean uplift pulses), or autocyclic processes (delta switching)
- Residence time in each mode affects ecosystem structure and mollusc community composition
- Mode 3 (decoupled) may represent storm events or short-lived perturbations that don't reset bottom water chemistry

---

## PART 7: SYNTHESIS & IMPLICATIONS (2-3 slides)

**NARRATIVE TRANSITION:**
> "Let's step back and ask: what does all of this mean? How do these geochemical patterns fit into our understanding of the Pebas mega-wetland as a whole?"

### Slide 16: Depositional Model
**Visual: Cross-section showing TAM and SC positions in Pebas system**

| Setting | Site | Characteristics | XRF Signature |
|---------|------|-----------------|---------------|
| **Lacustrine (distal)** | TAM | Stratified, dysoxic, Pachydon-dominated | Low Mn/Ti (0.18), high Ca/Ti (5.0) |
| **Fluvio-lacustrine (proximal)** | SC | Ventilated, variable, Diplodon + Pachydon | High Mn/Ti (0.37), high Ti (detrital) |

**TALKING POINTS:**
> "Imagine the Pebas mega-wetland as a series of interconnected environments: rivers bringing sediment and oxygen from the Andes, shallow embayments where mixing occurs, and deep, quiet lakes where stratification develops."

> "TAM sits in the distal, lacustrine part of this continuum. Water has traveled far from the sediment source - the Ti counts are lower because coarse detritus has already settled. Thermal stratification develops and persists."

> "SC is closer to the action. Higher Ti means more detrital input - we're nearer to fluvial channels. The system gets stirred up more often, breaking stratification and ventilating the bottom."

> "This isn't speculation. The XRF signatures confirm what paleontology suggested: a gradient from river-influenced to lake-dominated environments within the same formation, same time period, just 30 km apart."

**BACKGROUND - Pebas basin configuration:**
- Andean foreland basin: subsidence in the west, cratonic highlands in the east
- Rivers entered from both west (proto-Andes) and east (Guiana/Brazilian shields)
- Brackish marine incursions from the north (Caribbean/Venezuela) documented by Sr isotopes
- Basin evolved from open seaway to enclosed mega-wetland between ~23 and 8 Ma

### Slide 17: Drivers of Variability
**Proposed mechanisms for oxygenation events:**
1. **Andean runoff pulses** - Oxygenated river water entering lake
2. **Lake-level changes** - Shallowing increases mixing
3. **Orbital forcing** - Cyclic sediment deposition documented by Hoorn et al.

**Evidence**: Sr isotopes in mollusc shells identify Andean vs. cratonic water sources (Kaandorp 2005)

**TALKING POINTS:**
> "What drove these environmental swings? We can't point to a single cause, but three mechanisms are plausible, and they're not mutually exclusive."

> "Andean runoff pulses: The Andes were actively uplifting during the Miocene, generating sediment and runoff. Storm events or wet periods could have delivered cold, oxygen-rich water that plunged into the stratified lake and ventilated the bottom."

> "Lake-level changes: When lake level drops, you reduce the water column depth. A shallower lake is easier to mix - wind can reach the bottom. This would produce the ventilation events we see."

> "Orbital forcing: Hoorn showed cyclic sedimentation in Pebas cores at Milankovitch frequencies - precession, obliquity, eccentricity. These astronomical cycles modulate tropical rainfall and could drive periodic freshening or drying."

> "The Sr isotope evidence is compelling. Kaandorp showed that shells record shifts between Andean water (high ⁸⁷Sr/⁸⁶Sr from young volcanic rock) and cratonic water (low ⁸⁷Sr/⁸⁶Sr from old shields). The water sources changed through time."

**BACKGROUND - Testing these hypotheses:**
- Time series analysis of Mn/Ti could reveal periodicities consistent with orbital forcing
- Correlation with regional paleoclimate proxies (Amazon fan, Andes foreland) could test climate drivers
- Sequence stratigraphic analysis could test lake-level control
- Currently available data cannot definitively distinguish these mechanisms

### Slide 18: Conclusions
1. XRF Mn/Ti successfully detects dysoxic conditions independently documented by Pachydon ecology
2. Site-level contrast (TAM vs SC) matches lacustrine vs. fluvio-lacustrine interpretations from paleontology
3. Concordant Mn/Ti + Ca/Ti events provide strongest paleoenvironmental signals
4. Extreme transitions (833% Mn/Ti change) record discrete ventilation events in an otherwise stratified system

**TALKING POINTS:**
> "Four main conclusions, each building on the last."

> "First: XRF works for Pebas. Mn/Ti detects the same dysoxic conditions that Wesselingh identified from mollusc ecology 20 years ago. Independent validation is essential, and we have it."

> "Second: XRF reveals site differences that make paleontological sense. TAM and SC aren't just two points - they represent end members of an environmental gradient within the mega-wetland."

> "Third: When independent proxies agree, pay attention. Concordant Mn/Ti and Ca/Ti anomalies are your most robust signals. They resist the compositional closure critique because Ca and Mn have different source mechanisms."

> "Fourth: The system was dynamic. 833% transitions in Mn/Ti aren't noise - they're real ventilation events. The Pebas mega-wetland wasn't a uniform, static environment. It experienced dramatic environmental swings that shaped mollusc evolution."

> "Take-home message: XRF core scanning can extract paleoenvironmental information from Amazonian sediments that complements and validates paleontological interpretations."

**FINAL TALKING POINT (for Q&A preparation):**
> "What we still don't know: the absolute timing of events, the precise mechanisms driving variability, and how these local patterns connect to basin-wide Pebas evolution. These are the next questions."

---

## APPENDIX: Shiny App Settings for Figures

### Figure Generation Settings (http://127.0.0.1:6338)

**For site contrast figure:**
- Select: TAM-1-2-3B-A + SC-5-6-7ABC-B
- Proxies: Mn/Ti, Ca/Ti, Fe/Mn, Ti
- Sediment filter: OFF (show full section)
- Smoothing: 5
- Export: 16" × 10", 300 DPI

**For dramatic transitions:**
- Select: SC-1ABC-2-3C-A-RUN1 (show 350-400mm region)
- OR: TAM-5AB-6-7-B-RUN2 (show 630-670mm region)
- Smoothing: 3 (preserve sharp transitions)

**For concordant events:**
- Select: SC-3AB-4ABCD-RUN2-F (show 1240-1280mm region)
- Proxies: Mn/Ti, Ca/Ti only
- Highlight that these are INDEPENDENT proxies

---

## Key Intervals Quick Reference

### Must-Show Intervals:

| Priority | Section | Depth (mm) | Signal | Use For |
|----------|---------|------------|--------|---------|
| 1 | SC-1ABC-2-3C-A-RUN1 | 367 | +833% Mn/Ti | Most dramatic transition |
| 2 | SC-3AB-4ABCD-RUN2-F | 1258 | Concordant | Independent proxy validation |
| 3 | TAM-1-2-3B-A | 243-444 | Dysoxic baseline | Pachydon habitat validation |
| 4 | SC-5-6-7ABC-B | 704 | Ca/Ti = 48 | Extreme shell horizon |
| 5 | TAM-5AB-6-7-B-RUN2 | 264-285 | Sustained concordant | Best TAM oxic-carbonate |

### Best Sections for Each Site:

**TAM (lacustrine, dysoxic):**
1. TAM-1-2-3B-A (CV=1.28) - Most dynamic, best dysoxic intervals
2. TAM-5AB-6-7-B-RUN2 (CV=1.04) - Best transitions, concordant events
3. TAM-3A-4-5CDE-B (n=191) - Most measurements, clear carbonate pulses

**SC (fluvio-lacustrine, variable):**
1. SC-1ABC-2-3C-A-RUN1 (CV=1.24) - Most dramatic transitions
2. SC-5-6-7ABC-B (n=170) - Most measurements, extreme carbonate
3. SC-3AB-4ABCD-RUN2-F - Best concordant events

---

## COPY-PASTE TABLES FOR SLIDES

### Table 1: Pebas System Overview
```
PEBAS MEGA-WETLAND
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Area            >1 million km²
Duration        ~23–8 Ma (15 million years)
Age of cores    ~12 Ma (Middle Miocene)
Mollusc zones   MZ7–MZ8 (peak diversity)
Endemic species 158 molluscs (72% endemic)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
```

### Table 2: ITRAX Scanner Settings
```
XRF CORE SCANNING PARAMETERS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
X-ray tube      Mo (molybdenum)
Voltage         30 kV
Current         45 mA
Step size       2 mm
Integration     10 sec/point
Detector        Si-drift (SDD)
Output          Counts per second (cps)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
```

### Table 3: Element Behavior
```
ELEMENT GEOCHEMISTRY
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Element  Source      Redox?   Role
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Mn       Authigenic  YES      Primary O₂ tracer
Fe       Mixed       Partial  Supporting indicator
Ca       Biogenic    NO       Shell/productivity
Ti       Detrital    NO       Normalizer (inert)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
```

### Table 4: Proxy Interpretation
```
XRF PROXY FRAMEWORK
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Proxy    High Value           Low Value          Independent?
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Mn/Ti    Oxic (O₂ present)    Dysoxic (no O₂)    ✓ PRIMARY
Ca/Ti    Shell-rich           Shell-poor         ✓ YES
Fe/Mn    Reducing (anoxic)    Less reducing      ✗ shares Mn
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
```

### Table 5: Site Comparison
```
TAM vs SC: SITE CONTRAST
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Parameter        TAM       SC        Ratio    Meaning
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Mn/Ti median     0.18      0.37      2.1×     SC more oxic
Ca/Ti median     5.04      2.94      1.7×     TAM more shells
Ti (cps)         9,510     13,587    1.4×     SC more detrital
Mn (cps)         1,531     4,804     3.1×     SC more Mn precip
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
n measurements   876       1,062
Sections         11        16
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
```

### Table 6: Site Environments
```
DEPOSITIONAL SETTINGS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Site  Setting              Characteristics         XRF Signature
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
TAM   Lacustrine (distal)  Stratified, dysoxic     Low Mn/Ti, high Ca/Ti
                           Pachydon-dominated
SC    Fluvio-lacustrine    Ventilated, variable    High Mn/Ti, high Ti
      (proximal)           Mixed mollusc fauna
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
```

### Table 7: Key Intervals
```
HIGHLIGHT INTERVALS FOR PRESENTATION
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Section              Depth    Signal           Story
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
SC-1ABC-2-3C-A-RUN1  367 mm   +833% Mn/Ti      Most dramatic O₂ event
SC-3AB-4ABCD-RUN2-F  1258 mm  +5.8σ Mn/Ti      Best concordant signal
TAM-1-2-3B-A         243-444  Mn/Ti = 0.078    Extreme dysoxia
SC-5-6-7ABC-B        704 mm   Ca/Ti = 48       Massive shell bed
TAM-5AB-6-7-B-RUN2   264 mm   Concordant       Best TAM oxic event
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
```

### Table 8: Environmental Modes
```
THREE ENVIRONMENTAL STATES
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Mode  Mn/Ti  Ca/Ti  Environment          Best Example
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
1     HIGH   HIGH   Oxic + productive    SC-3AB-4ABCD-RUN2-F @ 1258
2     LOW    varies Dysoxic baseline     TAM-1-2-3B-A @ 243-444
3     LOW    HIGH   Shells + reducing    TAM-1-2-3B-B @ 1204
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
```

### Table 9: Isotope Validation
```
BIVALVE ISOTOPE COMPARISON (Kaandorp 2005)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Bivalve    Habitat      δ¹⁸O Pattern         XRF Match
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Pachydon   Lacustrine   Low amplitude        Stable low Mn/Ti (TAM)
           (dysoxic)    Irregular cycles     Stratified bottom water
Diplodon   Fluvial      High amplitude       Variable Mn/Ti (SC)
           (oxic)       Seasonal cycles      Ventilated seasonally
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
```

### Table 10: Mn Redox Chemistry
```
MANGANESE REDOX BEHAVIOR
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Condition    Mn State    Solubility    Sediment Effect
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
OXIC         Mn⁴⁺        Insoluble     Precipitates as MnO₂
(O₂ > 0.2)   (oxidized)  (trapped)     → HIGH Mn/Ti

DYSOXIC      Mn²⁺        Soluble       Diffuses into water
(O₂ < 0.2)   (reduced)   (mobile)      → LOW Mn/Ti
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Threshold: ~200 mV (Eh), ~0.2 ml/L O₂
```

### Table 11: Conclusions Summary
```
KEY FINDINGS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#   Finding                                    Evidence
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
1   XRF detects dysoxia                        Pachydon + low Mn/Ti match
2   Sites = different O₂ regimes               TAM 0.18 vs SC 0.37 Mn/Ti
3   Concordant proxies = robust signal         Ca/Ti + Mn/Ti independent
4   System was dynamic                         833% Mn/Ti transitions
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
```

### Table 12: Quick Numbers
```
NUMBERS TO REMEMBER
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
~12 Ma          Age of cores
~2,000          XRF measurements
2 mm            Spatial resolution
~3-5 yr         Temporal resolution
0.18 vs 0.37    TAM vs SC Mn/Ti
2.1×            SC/TAM oxygen difference
833%            Largest Mn/Ti transition
+5.8σ           Strongest anomaly
30 km           Distance between sites
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
```

### Table 13: Drivers of Variability
```
PROPOSED MECHANISMS FOR O₂ EVENTS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Mechanism          Process                  Evidence
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Andean runoff      Cold O₂-rich water       Sr isotopes in shells
                   plunges into lake
Lake-level drop    Shallowing enables       Sequence stratigraphy
                   wind mixing
Orbital forcing    Precession cycles        Cyclic sedimentation
                   modulate rainfall        (Hoorn et al.)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
```

### Table 14: Lithology Color Guide
```
READING CORE COLORS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Color           Chemistry       Environment
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Dark gray/black High organic    Reducing (dysoxic)
                Low Mn/Ti
Light gray/tan  Low organic     Oxidized (oxic)
                High Mn/Ti
White bands     High Ca         Shell-rich horizons
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
```

---

## References for Slides

- Aitchison, J. (1986). The Statistical Analysis of Compositional Data. Chapman and Hall.
- Hoorn, C. et al. (2010). The development of the Amazonian mega-wetland. In: Amazonia: Landscape and Species Evolution.
- Kaandorp, R.J.G. et al. (2005). Ecological implications from geochemical records of Miocene Western Amazonian bivalves. J. South American Earth Sciences.
- Wesselingh, F.P. et al. (2002). Lake Pebas: a palaeoecological reconstruction. Cainozoic Research.
- Wesselingh, F.P. (2006a). Evolutionary ecology of the Pachydontinae. Scripta Geologica.
