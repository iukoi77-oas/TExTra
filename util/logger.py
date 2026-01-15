import time

def get_args_info():
    print()

def get_current_time():
    formatted_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    return formatted_time

def log_message(level, message, bold=False, color="default"):
    colors = {
        "default": "\033[0m",
        "success": "\033[1;32m",   # Green
        "step": "\033[1;36m",   # Cyan
        "warning": "\033[1;33m",# Yellow
        "error": "\033[1;31m" # Red
    }
    start_color = colors.get(color, colors["default"])
    end_color = "\033[0m"
    bold_code = "\033[1m" if bold else ""
    
    print(f"{start_color}{level}{end_color} {bold_code}{message}{end_color}")
