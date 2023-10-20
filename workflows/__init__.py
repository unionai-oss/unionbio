import yaml
import logging
from pathlib import Path

# Determine project root
ROOT_PATH = Path(__file__).resolve().parent

# Import config
CONFIG_PATH = ROOT_PATH.joinpath('config.yaml')
with open(CONFIG_PATH, 'r') as f:
    config = yaml.load(f, Loader=yaml.FullLoader)

# Setup the logger
logger = logging.getLogger(__name__)
console_handler = logging.StreamHandler()
console_handler.setFormatter(logging.Formatter("[%(asctime)s %(levelname)s %(name)s] %(message)s"))
logger.addHandler(console_handler)
logger.setLevel(logging.DEBUG)