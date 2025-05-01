# logging_config.py
import logging

# Configure the logger
def setup_logger():
    logger = logging.getLogger("app_logger")  # Use a consistent name for the logger
    logger.setLevel(logging.DEBUG)

    # Create handler
    handler = logging.StreamHandler()  # StreamHandler logs to stdout
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)

    # Avoid adding multiple handlers if the logger is reused
    if not logger.handlers:
        logger.addHandler(handler)

    return logger

# Call setup_logger at the module level to ensure configuration
setup_logger()