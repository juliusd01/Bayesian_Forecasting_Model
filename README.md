# Election Forecast for Bundestagswahl 2025 using Bayesian Statistics

Main results for the national forecast:

![Vote Share Forecast](results/5_poll_modell/posterior_5poll_model.png)
![Vote Share Forecast](results/5_poll_modell/posterior_vote_share_5poll.png)
![Vote Share Forecast](results/voter_preferences.png)

Forecasts for each of the 16 states can be found in the results folder.


## Package overview
```txt
📦Term_paper
 ┣ 📂Resources_all                # Data sources used
 ┣ 📂results                      # The folder where all results are saved when running the scripts
 ┣ 📜2step.py                     # Script to generate main results 
 ┣ 📜2step_model.ipynb            # Notebook where it's easier to follow single steps, includes additional plots
 ┣ 📜5poll_model.ipynb            # Notebook for running the 5 poll model once
 ┣ 📜5poll_model.py               # Script to generate graphical and numerical results for the 5 poll model
 ┣ 📜README.md
 ┣ 📜beta_sensitivity.py          # Testing sensitivity of beta in the Gamma distribution in nations prior
 ┣ 📜c_sensitivity.py             # Testing the scaling factor in the Dirichlet distribution
 ┣ 📜create_legend_for_parties.py
 ┣ 📜environment.yml
 ┣ 📜license
 ┣ 📜plot_voter_preferences.py    # Script to generate Figure 11
 ┗ 📜requirements.txt
 ```

## Replicate results

- Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/products/distribution)

### Setup Environment
1. Clone or download this repository
2. Navigate to the project directory:
   ```bash
   cd /path/to/Bayesian_Forecasting_Model
   ```

3. Create the conda environment from the `environment.yml` file:
   ```bash
   conda env create -f environment.yml
   ```

4. Activate the environment:
   ```bash
   conda activate bayesian-forecast
   ```

### Run the Analysis
Execute the main scripts to generate results:
```bash
python 2step.py             
python 5poll_model.py 
python plot_voter_preferences.py
```

Results will be saved in the `results/` folder.






