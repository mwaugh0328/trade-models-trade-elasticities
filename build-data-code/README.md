Readme file to describe data construction. Still in construction. 

The main driver file is the ``make_price_data.jl`` which brings in the price and trade flow data for years 2004, 2011 and outputs ``pricegap-df-year.csv`` and ``tradeshare-df-year.csv`` which are used in the estimation. For the year 2017, the 2017 files are created with ``make2017_data.jl`` file.

Below are more details on each different element. 

- **Gross Output Data** This is simple, the file ``clean-grossoutput-data.ipynb'' grabs the UNIDO data for the three years 2004, 2011, 2017 and organizes to construct the trade share matrix.

- **ICP price data** TODO. In the folders ``Year-price-gap`` folder are the files for the ICP price data and then the file ``make_price_data.jl`` adjusts the price gap files for use later.

- **Trade data** For the year 2004, 2011 the trade flow data is from Feenstra's world trade flows (these do not appear to be publicly available anymore, will try with baci in the future). These are modified in the ``make_price_data.jl``or ``make_tradeshare_data.jl`` file. For the year 2017, the trade flow data is BACI and this is cleaned / shapped using ``clean-trade-data.ipynb`` python notebook.



