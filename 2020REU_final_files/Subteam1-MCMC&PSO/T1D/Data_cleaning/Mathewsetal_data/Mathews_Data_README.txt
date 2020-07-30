This data was sent to the DePillis+Shtylla 2019 REU team from Mathews-2015. The citation for the paper this data was used in is in APA format is here:

Mathews, C. E., Xue, S., Posgai, A., Lightfoot, Y. L., Li, X., Lin, A., … Atkinson, M. A. (2015). Acute Versus Progressive Onset of Diabetes in NOD Mice: Potential Implications for Therapeutic Interventions in Type 1 Diabetes. Diabetes, 64(11), 3885–3890. [https://doi.org/10.2337/db15-0449](https://doi.org/10.2337/db15-0449)

---

This data is time-series data tracking female NOD mice. Blood glucose readings were taken every two days, 
unless in the case that the blood sugar exceeded 250 mg/dL, in which case the blood glucose level was taken 24 hours later. 
Diabetes was diagnosed when two successive readings of blood sugar were above 250 mg/dL. 

Glucose levels are reported starting 24 days before diagnosis of diabetes. Glucose started being taken at 8 weeks 
(not necessarily reported here). 

---

Treatment arms: when diabetes was diagnosed, the mice was implanted with a subcutaneous insulin pellet. Data was taken for 
the next 120 days after the insulin pellet was inserted, or until two readings above 400 mg/dL were seen. Disease reversal 
was defined to be reaching a blood glucose level of less than 250 mg/dL after 120 days of the insulin treatment being started. 

---

Mice types: after data was taken, mice were classified according to the characteristics of their onset. If the mouse had 
steady euglycemia (glucose levels in the healthy range) before having a sudden onset of hyperglycemia, the mouse was 
categorized as having an acute form of the disease. Conversely, if the mouse simply had an 'excursion' above 200 mg/dL 
(meaning there was a reading of 200 or above, and at least one subsequent reading before diagnosis with diabetes was below this number), 
the diagnosis of the disease was categorized to be progressive. 

---

In XY time series tables (Progressive individual animals shown; acute+progressive; in_anm_acute; acute), each row is in an 
individual subject. The columns represent the day that the measurement was taken. Group A corresponds with -24, meaning that this 
glucose reading was taken 24 days before the diagnosis of diabetes; Group B corresponds with -23, meaning that this glucose reading was 
taken 23 days before the diagnosis of diabetes; etc. 

---

In mean / SD / N tables (acute+progresive_seperated; LinReg.Acute+Prog), each row is defined as the aggregate glucose 
levels for the date before diagnosis with diabetes. The mean and standard deviation, as well as the number of measurements 
taken for that day (n) are reported. 

---

A few rows in the time series data have an asterisk that do not let typical spreadsheet analytics software to process the data; 
we assume here that this data must be nullified for some reason (this is to be verified). 

---

The 2019 REU team working on this project were supported partially or fully by NSF Award 1757952. Should you have any questions about 
the project in totality, please contact Blerta Shtylla (Pomona College) or Lisette de Pillis (Harvey Mudd College).

Should you have questions for the author of this file, please contact Matt Matusiewicz by email at matt.matu@berkeley.edu or at 949.870.8373.


Notes from Shtylla converstation with A. Posgai on 08/12/2019:


The data that we have is what they also have on their end. Amanda has some original records that contain the actual day of onset for her studies. One issue is that this paper reports mice from a lot of other studies where they only kept records sporadically (it was not done by Amanda). She however noted that our graphs seem a lot more sparse than what their data looks like on her prism file. She said that they ploted the graphs using Prism. She also noted that we have a one day shift in our data files (so check our file transfer for this). Blerta will look at the plotting another time to make sure that there are no weird issues related to how the data is read into matlab plotting. Another thing to do is to try plotting using Prism.
In addition, Amanda said that she might have the onset day for her mice recorded so we can match at least some of the mice using the ID in the spreadsheets. She said that her protocol was what Clay sent to us, namely:

Blood Glucose Measures Prior to Onset
Nonfasting morning blood glucose concentrations were measured by tail vein prick using a OneTouch Ultra 2 blood glucose meter (LifeScan Inc.) every other day starting at 8 weeks of age. Any reading of blood glucose >250 mg/dL was followed by a test 24 h later. Diabetes was diagnosed by two successive readings of blood glucose >250 mg/dL (19). The blood glucose levels are reported starting at 24 days prior to diagnosis, a time period considered practical for evaluation based on the starting age of study (i.e., 8 weeks) relative to the presumed time period of when the earliest cases of diabetes would occur.













