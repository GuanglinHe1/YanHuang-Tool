#!/usr/bin/env bash
load_config() {
    local config_text
    config_text=$(Rscript "$PROJECT_ROOT/script/load_config.R" "$1") || return
    eval "$config_text"
}
