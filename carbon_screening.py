import cobra
import pandas as pd
def carbon_source_screening_rxn(model, target, default_carbon='EX_glc__D_e',target_rxn=None):
    """
    Screen the model for target production with different carbon sources and return the results.
    
    Parameters:
    - model: The metabolic model to screen.
    - target: The target metabolite to maximize (e.g., 'GALUi' for UDP-glucose).
    - target_rxn: The reaction ID for the target metabolite (e.g., 'GALUi').
    
    Returns:
    - A DataFrame with the results of the screening.
    """
    target_flux=pd.DataFrame(columns=['Exchange Reaction','metabolite','GrowthRate','Target Flux','Target Max Flux'])

    # QUitamos la fuente de carbono para ir rotando con las diferentes fuentes de carbono
    model.reactions.get_by_id(default_carbon).lower_bound = 0.0
    for rxn in model.exchanges:
        model.objective="BIOMASS_KT2440_WT3"
        old_lower=model.reactions.get_by_id(rxn.id).lower_bound
        for metabolite in rxn.metabolites:
            # Check if the carbon is in the metabolite formula
            if 'C' in metabolite.elements.keys():
                print(f"Testing exchange reaction: {rxn.id}")
                # Set the exchange reaction to be active
                # Normalizados el uptake según el número de carbonos que tenga
                carbon_count = sum(metabolite.elements.get('C', 0) for metabolite in rxn.metabolites)
                carbon_uptake= -36.0 / carbon_count if carbon_count > 0 else -6.0 # Normalizamos la fuente de carbono OJO!!!!!!!!
                print(f"Carbon uptake for {rxn.id}: {carbon_uptake}")
                # Sacamos el flujo de UDP-glucosa sin maximizar
                model.reactions.get_by_id(rxn.id).lower_bound = carbon_uptake
                solution= model.optimize()
                growth_rate = solution.objective_value
                #Flujo de UDP-glucosa sin maximizar
                # testear si el modelo crece y su solucion no es infeasible
                if solution.status == 'infeasible':
                    print(f"Model is infeasible with {rxn.id} and {metabolite.name}. Skipping...")
                    continue
                # Cogemos el flujo hacia el metabolito objetivo sin maximizar [solo la primera reacción que lo produce que es el maximo]
                target_flux_tmp= float(model.metabolites.get_by_id(target).summary().producing_flux.flux.iloc[0])
                # si hay target_rxn, devuelve el flujo de esa reacción sin maximizar
                if target_rxn is not None:
                    target_flux_tmp= float(solution.fluxes.get(target_rxn, 0.0))
                # Set the objective to maximize the target metabolite production
                if target_rxn is not None:
                    model.objective = target_rxn
                else: # If no target reaction is specified, use the exchange reaction as the target
                    # Obtener las reacciones asociadas al metabolito objetivo y devolver la que sea una exchange reaction
                    if target in model.metabolites:
                        target_reactions = model.metabolites.get_by_id(target).reactions
                        for reac in target_reactions:
                            if reac.id.startswith('EX_'):
                                model.objective = reac.id
                                break
                sol = model.slim_optimize()

                # Print the results
                #sol=float(model.metabolites.get_by_id(target).summary().producing_flux.flux.iloc[0])
            print(f"Exchange Reaction: {rxn.id}, metabolite:{metabolite.name},GrowthRate: {growth_rate}, Target Flux: {target_flux_tmp}, Target Max Flux: {sol}")
            # Create a DataFrame to store the results
            tmp_df= pd.DataFrame({
                'Exchange Reaction': [rxn.id],
                'metabolite': [metabolite.name],
                'GrowthRate': [growth_rate],
                'Reaction Flux': [target_flux_tmp],
                'Reaction Max Flux': [sol]
            })
            # Append the results to the udp_flux DataFrame
            
            target_flux = pd.concat([target_flux, tmp_df], ignore_index=True)

            # Reset the exchange reaction bounds for the next iteration
            model.reactions.get_by_id(rxn.id).lower_bound = old_lower     