import yaml


class DotDict(dict):
    def __getattribute__(self, name):
        return self[name]


def load_config(path="configs/config.yaml"):
    with open(path, "r") as config_file:
        config = yaml.safe_load(config_file)

    for key, value in config.items():
        if isinstance(value, dict):
            config[key] = DotDict(value)

    return config
