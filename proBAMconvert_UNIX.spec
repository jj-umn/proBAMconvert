# -*- mode: python -*-

block_cipher = None


a = Analysis(['proBAM.py'],
             pathex=['/home/vladie/PycharmProjects/proBAMconvert'],
             binaries=[],
             datas=[('proBAMconvert_logo.gif', '.')],
             hiddenimports=[],
             hookspath=[],
             runtime_hooks=[],
             excludes=['sqlite3', 'sqlalchemy', 'matplotlib', 'scipy', 'IPython', 'PyQt4'],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          name='proBAMconvert_UNIX',
          debug=False,
          strip=False,
          upx=True,
          console=True )
