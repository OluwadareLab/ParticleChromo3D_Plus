import logging

def setup_logger(log_level = logging.DEBUG):
    logger = logging.getLogger('ParticleChromo3D')  # Create a named logger
    logger.setLevel(log_level)  # You can set this to DEBUG, WARNING, etc.
    
    # Create formatters
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    # Create file handler for logging to a file
    file_handler = logging.FileHandler('ParticleChromo3d.log')
    file_handler.setFormatter(formatter)

    # Create stream handler for logging to the console
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)

    # Add handlers to the logger
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)

    return logger