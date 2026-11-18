// Copyright (C) 2026 CRV Science Network (https://crvscience.com)
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://gnu.org>.

#![forbid(unsafe_code)]

use std::env;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

fn parse_version_quad(version: &str) -> (u16, u16, u16, u16) {
    let core = version.split('-').next().unwrap_or("0.1.0");
    let mut parts = core.split('.').map(|p| p.parse::<u16>().unwrap_or(0));
    let major = parts.next().unwrap_or(0);
    let minor = parts.next().unwrap_or(0);
    let patch = parts.next().unwrap_or(0);
    (major, minor, patch, 0)
}

fn find_windres() -> Option<PathBuf> {
    for candidate in ["x86_64-w64-mingw32-windres", "windres"] {
        if Command::new(candidate)
            .arg("--version")
            .output()
            .is_ok_and(|o| o.status.success())
        {
            return Some(PathBuf::from(candidate));
        }
    }
    None
}

fn find_rc_exe() -> Option<PathBuf> {
    if Command::new("rc.exe")
        .arg("/?")
        .output()
        .is_ok_and(|o| o.status.success() || !o.stdout.is_empty() || !o.stderr.is_empty())
    {
        return Some(PathBuf::from("rc.exe"));
    }

    let kits_root = Path::new(r"C:\Program Files (x86)\Windows Kits\10\bin");
    if let Ok(entries) = fs::read_dir(kits_root) {
        let mut versions: Vec<PathBuf> = entries
            .filter_map(|e| e.ok().map(|x| x.path()))
            .filter(|p| p.is_dir())
            .collect();
        versions.sort();
        versions.reverse();
        for ver_dir in versions {
            for arch in ["x64", "x86", "arm64"] {
                let cand = ver_dir.join(arch).join("rc.exe");
                if cand.is_file() {
                    return Some(cand);
                }
            }
        }
    }
    None
}

fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-changed=assets/causal-atlas.ico");

    if env::var("CARGO_CFG_TARGET_OS").as_deref() != Ok("windows") {
        return;
    }

    let manifest_dir = PathBuf::from(env::var("CARGO_MANIFEST_DIR").expect("CARGO_MANIFEST_DIR"));
    let out_dir = PathBuf::from(env::var("OUT_DIR").expect("OUT_DIR"));
    let ico_path = manifest_dir.join("assets").join("causal-atlas.ico");
    if !ico_path.is_file() {
        panic!("Missing required Windows icon asset: {}", ico_path.display());
    }

    let pkg_version = env::var("CARGO_PKG_VERSION").unwrap_or_else(|_| "0.2.0".to_string());
    let (v1, v2, v3, v4) = parse_version_quad(&pkg_version);
    let ico_rc_path = ico_path.to_string_lossy().replace('\\', "/");

    let rc_content = format!(
        r#"1 ICON "{ico_rc_path}"

1 VERSIONINFO
FILEVERSION {v1},{v2},{v3},{v4}
PRODUCTVERSION {v1},{v2},{v3},{v4}
FILEFLAGSMASK 0x3fL
FILEFLAGS 0x0L
FILEOS 0x40004L
FILETYPE 0x1L
FILESUBTYPE 0x0L
BEGIN
  BLOCK "StringFileInfo"
  BEGIN
    BLOCK "040904b0"
    BEGIN
      VALUE "CompanyName", "CRV Science Network"
      VALUE "FileDescription", "Causal Atlas - Pathway, Upstream Regulator & Causal Network Workbench"
      VALUE "FileVersion", "{pkg_version}"
      VALUE "InternalName", "causal-atlas"
      VALUE "LegalCopyright", "Copyright CRV Science Network (GPL-3.0-or-later)"
      VALUE "OriginalFilename", "causal-atlas.exe"
      VALUE "ProductName", "Causal Atlas"
      VALUE "ProductVersion", "{pkg_version}"
    END
  END
  BLOCK "VarFileInfo"
  BEGIN
    VALUE "Translation", 0x409, 1200
  END
END
"#
    );

    let rc_path = out_dir.join("app.rc");
    fs::write(&rc_path, rc_content).expect("write app.rc");

    let target_env = env::var("CARGO_CFG_TARGET_ENV").unwrap_or_default();
    if target_env == "gnu" {
        if let Some(windres) = find_windres() {
            let res_obj = out_dir.join("app_res.o");
            let status = Command::new(&windres)
                .arg("-i")
                .arg(&rc_path)
                .arg("-O")
                .arg("coff")
                .arg("-o")
                .arg(&res_obj)
                .status()
                .expect("invoke windres");
            if !status.success() {
                panic!("windres failed with status {status}");
            }
            println!("cargo:rustc-link-arg-bins={}", res_obj.display());
        } else {
            println!("cargo:warning=windres not found; skipping Windows PE resource embedding");
        }
    } else if target_env == "msvc" {
        if let Some(rc_exe) = find_rc_exe() {
            let res_file = out_dir.join("app.res");
            let status = Command::new(&rc_exe)
                .arg("/nologo")
                .arg("/fo")
                .arg(&res_file)
                .arg(&rc_path)
                .status()
                .expect("invoke rc.exe");
            if !status.success() {
                panic!("rc.exe failed with status {status}");
            }
            println!("cargo:rustc-link-arg-bins={}", res_file.display());
        } else {
            println!("cargo:warning=rc.exe not found; skipping Windows PE resource embedding");
        }
    }
}
