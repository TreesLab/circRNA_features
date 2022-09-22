import csv



def read_sample_list(sample_list_file):
    with open(sample_list_file, newline='', encoding='utf-8') as f_in:
        sample_list_reader = csv.DictReader(f_in, delimiter='\t')
        return list(sample_list_reader)

SAMPLES_LIST = read_sample_list(config['samples'])


def outputs_with_sample_info(output_templ, samples_list, keys_to_be_filled):
    outputs = [
        output_templ.format(*[sample[key] for key in keys_to_be_filled])
        for sample in samples_list
    ]
    return outputs


def get_all_sample_groups():
    sample_groups = list(
        set(
            (
                sample['dataset_id'],
                sample['species'],
                sample['tissue']
            )
            for sample in SAMPLES_LIST
        )
    )
    sample_groups = sorted(sample_groups, key=lambda sample: [sample[0], sample[2]])
    return sample_groups

ALL_SAMPLE_GROUPS = get_all_sample_groups()


def get_treatment_names(wildcard):
    samples = list(
        filter(
            lambda sample: (sample['dataset_id'] == wildcard.dataset_id) and \
                (sample['species'] == wildcard.species) and \
                (sample['tissue'] == wildcard.tissue),
            SAMPLES_LIST
        )
    )
    
    sample_order = ('totalRNA', 'RNaseR', 'polyA_minus')
    sorted_samples = sorted(samples, key=lambda sample: sample_order.index(sample['treatments']))

    treatment_names = [sample['treatments'] for sample in sorted_samples]

    return treatment_names


def get_sample_id_mapper(sample_id_mapping_file):
    sample_id_mapper = {}
    with open(sample_id_mapping_file, newline='', encoding='utf-8') as f_in:
        sample_id_reader = csv.DictReader(f_in, delimiter='\t')
        for sample in sample_id_reader:
            sample_id_mapper[sample['old_id']] = [sample['new_id'], int(sample['order'])]

        return sample_id_mapper

SAMPLE_ID_MAPPER = get_sample_id_mapper(config['sample_id_mapping'])
