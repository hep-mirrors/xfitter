import yaml

# Here you should place your YAML file
with open('data_incljet_antiktr04_incljetpt_fullrun2_v17_eta9.yaml', 'r') as input_file:
    yaml_data = yaml.load(input_file, Loader=yaml.FullLoader)

# Extract relevant information
independent_variables = yaml_data['independent_variables'][0]['values']
dependent_variables = yaml_data['dependent_variables'][0]['values']

# Open your output file
with open('output9.txt', 'w') as output_file:
    for i in range(len(independent_variables)):
        low = independent_variables[i]['low']
        high = independent_variables[i]['high']
        value = dependent_variables[i]['value']
        errors = dependent_variables[i]['errors']
        
        # Write output values 
        output_file.write(f'{low} {high} {value}')
        
        # Add errors
        for error in errors:
            plus = error['asymerror']['plus']
            minus = error['asymerror']['minus']
            output_file.write(f' {plus} {minus}')
        
        output_file.write('\n')
