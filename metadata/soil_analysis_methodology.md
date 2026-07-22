# Reference soil analysis methodology (ZALF Central Laboratory, Müncheberg)

The wet-chemistry reference values used to calibrate and validate the
spectral models (`db_soil_analysis_all.csv`) were determined at the
**ZALF Zentrallabor (Central Laboratory), Müncheberg, Germany**, following
their standard operating procedures. The methods relevant to this study are
summarised below.

## pH (H₂O)

Soil pH was measured potentiometrically in a soil–water suspension (standard
ZALF soil:solution ratio) using a calibrated pH electrode. The laboratory
also offers KCl and CaCl₂ suspensions as alternative extractants, but the
water-extractable value was used as the reference for this study, as it is
the most common reference for comparing spectrally predicted pH in soil
spectroscopy literature.

## Total Carbon (Ct) and Total Nitrogen (Nt)

Total carbon and total nitrogen were determined by **dry combustion
(elemental analysis)**: a dried, ground soil subsample is combusted at high
temperature in an elemental analyser, converting all carbon and nitrogen to
CO₂ and N₂/NOₓ, which are then quantified by thermal conductivity or
infrared detection. This method captures total C (organic + any inorganic
carbonate C present) and total N, and is the standard reference method
against which NIR-predicted C and N are benchmarked.

## Plant-available Phosphorus — Double Lactate (DL) extraction

Plant-available P was determined using the **Double Lactate (DL) method**
(calcium-lactate/acetic acid extraction, buffered to pH ≈ 3.6), a standard
German agronomic soil-testing procedure for estimating the readily
plant-available P fraction. Extracted P is quantified colorimetrically
(molybdate-blue method) or by ICP-OES, depending on the lab's routine
workflow. The **CAL method** (calcium-acetate-lactate, a closely related
buffered extraction) is offered by the laboratory as an alternative and was
not used here.

## Notes on comparability

- All reference analyses were carried out on the air-dried, sieved (< 2 mm)
  soil fraction, consistent with standard practice for pairing wet-chemistry
  values with diffuse reflectance spectra.
- Total C from dry combustion includes any carbonate-C present; where
  carbonates were negligible, Total C was treated as equivalent to soil
  organic carbon for modelling purposes.
- These are standard, widely used German/European agronomic and analytical
  methods (not proprietary to ZALF); the laboratory follows its own internal
  SOPs for sample preparation, but the underlying chemistry (dry combustion,
  DL extraction, potentiometric pH) is consistent with common soil-testing
  practice.
