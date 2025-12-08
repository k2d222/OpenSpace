{
  "additional_scripts": [
    "openspace.setPropertyValueSingle(\"Scene.Earth.Renderable.Layers.ColorLayers.noaa-sos-land-forests-map.Remove\", nil)",
    "openspace.setPropertyValueSingle(\"Scene.Earth.Renderable.Layers.ColorLayers.noaa-sos-land-forests-gain.Remove\", nil)",
    "openspace.setPropertyValueSingle(\"Scene.Earth.Renderable.Layers.ColorLayers.noaa-sos-land-forests-loss.Remove\", nil)",
    "openspace.setPropertyValueSingle(\"Scene.Earth.Renderable.Layers.ColorLayers.noaa-sos-land-forests.Remove\", nil)"
  ],
  "assets": [
    "base_blank",
    "scene/digitaluniverse/backgroundradiation",
    "scene/digitaluniverse/constellations",
    "scene/solarsystem/planets/earth/layers/colorlayers/blue_marble",
    "scene/solarsystem/planets/earth/layers/colorlayers/viirs_noaa21_temporal",
    "scene/solarsystem/planets/earth/layers/nightlayers/earth_at_night_2012",
    "scene/solarsystem/planets/earth/moon/layers/colorlayers/lola_clr_shade_sweden",
    "scene/solarsystem/planets/earth/moon/layers/colorlayers/wac_v1_sweden",
    "scene/solarsystem/planets/earth/noaa-sos/atmosphere/land_temp",
    "scene/solarsystem/planets/earth/noaa-sos/atmosphere/lightning",
    "scene/solarsystem/planets/earth/noaa-sos/land/agriculture-cropland",
    "scene/solarsystem/planets/earth/noaa-sos/land/flooding-major_floods",
    "scene/solarsystem/planets/earth/noaa-sos/land/forests",
    "scene/solarsystem/planets/earth/noaa-sos/land/top_quakes",
    "scene/solarsystem/planets/earth/noaa-sos/oceans/chlorophyll_model",
    "scene/solarsystem/planets/earth/noaa-sos/oceans/currents",
    "scene/solarsystem/planets/earth/noaa-sos/oceans/vector_winds",
    "scene/solarsystem/planets/earth/noaa-sos/overlays/capitals",
    "scene/solarsystem/planets/earth/noaa-sos/overlays/continent_names",
    "scene/solarsystem/planets/earth/noaa-sos/overlays/country_borders-white",
    "scene/solarsystem/planets/earth/noaa-sos/overlays/currents",
    "scene/solarsystem/planets/earth/noaa-sos/overlays/general_circulation",
    "scene/solarsystem/planets/earth/noaa-sos/overlays/ocean_names",
    "scene/solarsystem/planets/earth/noaa-sos/overlays/plate_boundary-color",
    "scene/solarsystem/planets/earth/noaa-sos/overlays/plate_names",
    "scene/solarsystem/planets/earth/noaa-sos/overlays/railroad",
    "scene/solarsystem/planets/earth/noaa-sos/overlays/rivers",
    "scene/solarsystem/planets/earth/noaa-sos/overlays/roads-white",
    "scene/solarsystem/planets/earth/noaa-sos/overlays/timezones",
    "scene/solarsystem/planets/jupiter/callisto/layers/colorlayers/callisto_texture",
    "scene/solarsystem/planets/jupiter/europa/layers/colorlayers/europa_texture",
    "scene/solarsystem/planets/jupiter/europa/layers/colorlayers/voyager_global_mosaic_sweden",
    "scene/solarsystem/planets/jupiter/ganymede/layers/colorlayers/ganymede_texture",
    "scene/solarsystem/planets/jupiter/io/layers/colorlayers/io_texture",
    "scene/solarsystem/planets/jupiter/layers/colorlayers/jupiter_texture",
    "scene/solarsystem/planets/mars/layers/colorlayers/mars_texture",
    "scene/solarsystem/planets/mars/layers/colorlayers/moc_wa_color_sweden",
    "scene/solarsystem/planets/mars/layers/colorlayers/mola_pseudo_color_sweden",
    "scene/solarsystem/planets/mercury/layers/colorlayers/messenger_mdis_sweden",
    "scene/solarsystem/planets/mercury/layers/colorlayers/messenger_mosaic2_sweden",
    "scene/solarsystem/planets/mercury/layers/colorlayers/messenger_shade_sweden",
    "scene/solarsystem/planets/neptune/layers/colorlayers/neptune_texture",
    "scene/solarsystem/planets/neptune/triton/layers/colorlayers/triton_voyager2_clrmosaic_globalfill_600m",
    "scene/solarsystem/planets/saturn/dione/layers/colorlayers/dione_texture",
    "scene/solarsystem/planets/saturn/enceladus/layers/colorlayers/enceladus_texture",
    "scene/solarsystem/planets/saturn/enceladus/layers/colorlayers/global_mosaic_100m_hpf_sweden",
    "scene/solarsystem/planets/saturn/iapetus/layers/colorlayers/iapetus_texture",
    "scene/solarsystem/planets/saturn/layers/colorlayers/saturn_texture",
    "scene/solarsystem/planets/saturn/rhea/layers/colorlayers/rhea_texture",
    "scene/solarsystem/planets/saturn/tethys/layers/colorlayers/tethys_texture",
    "scene/solarsystem/planets/saturn/titan/default_layers",
    "scene/solarsystem/planets/saturn/titan/layers/colorlayers/cassini_iss_global_mosaic_4km_sweden",
    "scene/solarsystem/planets/saturn/titan/layers/colorlayers/cassini_sar_hisar_global_mosaic_351m_sweden",
    "scene/solarsystem/planets/saturn/titan/titan",
    "scene/solarsystem/planets/uranus/layers/colorlayers/uranus_texture",
    "scene/solarsystem/planets/venus/layers/colorlayers/venus_texture",
    "scene/solarsystem/sun/layers/colorlayers/sun_texture"
  ],
  "camera": {
    "altitude": 17000000.0,
    "anchor": "Earth",
    "latitude": 58.5877,
    "longitude": 16.1924,
    "type": "goToGeo"
  },
  "delta_times": [
    1.0,
    5.0,
    30.0,
    60.0,
    300.0,
    1800.0,
    3600.0,
    43200.0,
    86400.0,
    604800.0,
    1209600.0,
    2592000.0,
    5184000.0,
    7776000.0,
    15552000.0,
    31536000.0,
    63072000.0,
    157680000.0,
    315360000.0,
    630720000.0
  ],
  "mark_nodes": [
    "Earth"
  ],
  "meta": {
    "author": "OpenSpace Team",
    "description": "Default OpenSpace Profile. Adds Earth satellites not contained in other profiles",
    "license": "MIT License",
    "name": "Default",
    "url": "https://www.openspaceproject.com",
    "version": "1.0"
  },
  "properties": [
    {
      "name": "Scene.*Label.Renderable.Enabled",
      "type": "setPropertyValue",
      "value": "false"
    },
    {
      "name": "Scene.*Trail.Renderable.Enabled",
      "type": "setPropertyValue",
      "value": "false"
    },
    {
      "name": "Scene.Earth.Renderable.Layers.Overlays.*.Enabled",
      "type": "setPropertyValue",
      "value": "false"
    },
    {
      "name": "Scene.Earth.Renderable.Layers.ColorLayers.*.Enabled",
      "type": "setPropertyValue",
      "value": "false"
    },
    {
      "name": "Scene.Earth.Renderable.Layers.NightLayers.*.Enabled",
      "type": "setPropertyValueSingle",
      "value": "false"
    },
    {
      "name": "Scene.Earth.Renderable.Layers.ColorLayers.Blue_Marble.Enabled",
      "type": "setPropertyValueSingle",
      "value": "true"
    },
    {
      "name": "Scene.Moon.Renderable.Layers.ColorLayers.Lola_Clr_Shade_Sweden.Enabled",
      "type": "setPropertyValueSingle",
      "value": "false"
    },
    {
      "name": "Scene.Mars.Renderable.Layers.ColorLayers.MOLA_Pseudo_Color_Sweden.Enabled",
      "type": "setPropertyValueSingle",
      "value": "false"
    },
    {
      "name": "Scene.Sun.Renderable.Enabled",
      "type": "setPropertyValueSingle",
      "value": "true"
    },
    {
      "name": "Scene.*.Renderable.PerformShading",
      "type": "setPropertyValue",
      "value": "false"
    },
    {
      "name": "Scene.Mercury.Renderable.Layers.ColorLayers.Messenger_SHADE_Sweden.Enabled",
      "type": "setPropertyValueSingle",
      "value": "false"
    },
    {
      "name": "ScreenSpace.*.Enabled",
      "type": "setPropertyValue",
      "value": "false"
    }
  ],
  "time": {
    "is_paused": false,
    "type": "relative",
    "value": "-1d"
  },
  "version": {
    "major": 1,
    "minor": 4
  }
}