import yaml

# Apro il file YAML in lettura
with open('data_incljet_antiktr04_incljetpt_fullrun2_v17_eta9.yaml', 'r') as input_file:
    yaml_data = yaml.load(input_file, Loader=yaml.FullLoader)

# Estrai le informazioni dai dati YAML
independent_variables = yaml_data['independent_variables'][0]['values']
dependent_variables = yaml_data['dependent_variables'][0]['values']

# Apro il file di output in scrittura
with open('output9.txt', 'w') as output_file:
    # Ciclo attraverso le coppie di valori independent_variables e dependent_variables
    for i in range(len(independent_variables)):
        low = independent_variables[i]['low']
        high = independent_variables[i]['high']
        value = dependent_variables[i]['value']
        errors = dependent_variables[i]['errors']
        
        # Scrivo i valori nel formato richiesto nel file di output
        output_file.write(f'{low} {high} {value}')
        
        # Aggiungo i valori di errore
        for error in errors:
            plus = error['asymerror']['plus']
            minus = error['asymerror']['minus']
            output_file.write(f' {plus} {minus}')
        
        # Vado a capo alla fine di ogni riga
        output_file.write('\n')
