
#--------------------------------------------------------------------------------
#climatologies
############################################################################################################
# Read the CSV file
# Name of the file in the current directory
file_name = "climatologies.csv"

# Read the CSV file
clima_df = CSV.read(file_name, DataFrame; delim=';', decimal=',')
dropmissing!(clima_df)

# Remove the 366th row
clima_df = clima_df[1:365, :]

# Repeat the DataFrame 50 times
repeated_df = vcat([clima_df for _ in 1:160]...)

clima_X_L = Vector(repeated_df[!, :JL_mean]) #16071 elements from 1.1.1975 to 31.12.2018
clima_X_tot = Vector(repeated_df[!, :Jtot_mean]) #16071 elements from 1.1.1975 to 31.12.2018
clima_temp = Vector(repeated_df[!, :thetao_mean])
clima_Mf0 = Vector(repeated_df[!, :mf0_pil_mean])
clima_Mf1 = Vector(repeated_df[!, :mf1_pil_mean])
clima_Mf2 = Vector(repeated_df[!, :mf2_pil_mean])
clima_Mf3 = Vector(repeated_df[!, :mf3_pil_mean])
clima_Mf4 = Vector(repeated_df[!, :mf4_pil_mean])