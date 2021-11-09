# The Effectiveness of Various Precautions for Pandemics

## Motivation
Test the effectiveness of various precautions for pandemics such as wearing masks, social distancing, quarantine, and so on. See how effective these preventive measures is. Find a best way to control the pandemic.

## Model
### SIR Model
There are three kinds of people in model.
- Susceptible (S): blue color
- Infectious (I): red color
- Recovered or Removed (died) (R): gray color\

Initialized parameters:
- `infection_radius` = 0.5
- `probability_of_infection` = 0.5

## Mechanism
- **Each individual walks in random.**
- **Individuals can move between different communities.**\
  It has 6 communities in the model.\
  Initialized parameters:
  - `travel_rate` = 0.3 (between communities)
  - `trigger_case_rate` = 0.05\
  Start decreasing the travel rate after the `trigger_case_rate`.
- **Infection Identification:**\
  We set 0.2 units of time in the program is one day in real world.\
  Initialized parameters:
  - `percentage_of_no_symptom` = 0.3
  - `infectious_period` = 4.2 (21 days)
- **Quarantine**\
  Initialized parameters:
  - `quarantine_mode` = `false`
  - `percentage_of_quarantine` = 0.7
- **Social distancing**\
  Initialized parameters:
  - `social_distanc_factor` = sqrt(0.5)\
  The value in the square root means the distance between people in real world.\
  When this factor is larger, the social distance will be closer.
  - `percentage_of_social_distancing` = 0.7

## Transmission Level Evaluation
- **Rate of spread:**\
  `R = (number of infectious people in current moment) / (number of infected people in current moment)`
  - R > 1: Infected people are exponentially increasing.
  - R = 1: The pandemic is under control.
  - R < 1: The pandemic is slowing down.

## Conclusion
The precedence of these precautions from highest to lowest:
1. **Wearing masks** and **quaratine**
2. **Larger social distance** and **percentage of people following social distancing**
3. **Reducing travel rate**
